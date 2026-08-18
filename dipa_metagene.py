#!/usr/bin/env python3
"""
Create exon-only metagene profiles around significant UP dIPAs.
"""

import argparse
import json
import os
import re
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import gffutils

# Restricted HPC nodes do not always allow applications to write under $HOME.
# Use node-local temporary directories unless the user already configured them.
temporary_directory = Path(tempfile.gettempdir())
matplotlib_config_directory = (
    temporary_directory / f"dipa_metagene_matplotlib_{os.getuid()}"
)
general_cache_directory = (
    temporary_directory / f"dipa_metagene_cache_{os.getuid()}"
)
matplotlib_config_directory.mkdir(parents=True, exist_ok=True)
general_cache_directory.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_config_directory))
os.environ.setdefault("XDG_CACHE_HOME", str(general_cache_directory))

import matplotlib

# HPC compute nodes usually do not have a graphical display.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig


REQUIRED_APA_COLUMNS = [
    "gene_symbol",
    "PASid",
    "RED",
    "p_adj",
    "APAreg",
    "APAreg_original",
    "chr",
    "start",
    "end",
    "browser_start_1based",
    "pas_coordinate_match",
]

REQUIRED_SAMPLE_COLUMNS = [
    "sample_id",
    "condition",
    "replicate",
    "role",
    "pair_id",
    "bigwig",
]

TRUE_TEXT_VALUES = {"true", "t", "1", "yes", "y"}

COLORBLIND_COLORS = [
    "#0173B2",
    "#DE8F05",
    "#029E73",
    "#D55E00",
    "#CC78BC",
    "#CA9161",
    "#FBAFE4",
    "#949494",
    "#ECE133",
    "#56B4E9",
]


def add_exclusion(excluded_rows, source_row, stage, reason, details=""):
    """Add one row to the exclusion QC table."""

    if hasattr(source_row, "to_dict"):
        source_values = source_row.to_dict()
    else:
        source_values = dict(source_row)

    clean_values = {
        key: value
        for key, value in source_values.items()
        if not str(key).startswith("_")
    }
    clean_values["exclusion_stage"] = stage
    clean_values["exclusion_reason"] = reason
    clean_values["details"] = details
    excluded_rows.append(clean_values)


def parse_pas_strand(pas_id):
    """Extract a legacy + or - label from a PASid when one is present."""

    match = re.search(r":([+-]):[^:]+$", str(pas_id).strip())
    if match is None:
        return None
    return match.group(1)


def normalize_gene_id(gene_id):
    """Remove a version suffix so ENSMUSG... and ENSMUSG....1 can match."""

    return str(gene_id).strip().split(".")[0]


def create_or_open_gtf_database(gtf_path, database_path, rebuild_database):
    """Create the reusable gffutils database, or open it when it already exists."""

    if rebuild_database and database_path.exists():
        database_path.unlink()

    if not database_path.exists():
        print(f"Creating GTF database: {database_path}")
        print("This can take several minutes for a full GENCODE annotation.")
        gffutils.create_db(
            str(gtf_path),
            dbfn=str(database_path),
            force=True,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True,
            verbose=False,
        )
    else:
        print(f"Reusing GTF database: {database_path}")

    return gffutils.FeatureDB(str(database_path), keep_order=True)


def prepare_cdna_model(apa_row, transcript, gtf_database, min_side_bp):
    """Build upstream and downstream exon lists for one selected dIPA."""

    transcript_exons = list(
        gtf_database.children(
            transcript,
            featuretype="exon",
            order_by="start",
        )
    )

    if len(transcript_exons) < 2:
        return None, "canonical_transcript_has_fewer_than_two_exons", ""

    exons = []
    for exon in transcript_exons:
        exon_id_values = exon.attributes.get("exon_id", [])
        exon_id = exon_id_values[0] if exon_id_values else exon.id
        exons.append(
            {
                "exon_id": exon_id,
                "start_1based": int(exon.start),
                "end_1based": int(exon.end),
                "start_0based": int(exon.start) - 1,
                "end_0based": int(exon.end),
            }
        )

    exons.sort(
        key=lambda exon: exon["start_1based"],
        reverse=(transcript.strand == "-"),
    )

    dipa_position = int(apa_row["_browser_start_numeric"])
    transcript_start = min(exon["start_1based"] for exon in exons)
    transcript_end = max(exon["end_1based"] for exon in exons)

    if not transcript_start <= dipa_position <= transcript_end:
        details = (
            f"dIPA={dipa_position}; transcript_span="
            f"{transcript_start}-{transcript_end}"
        )
        return None, "dipa_outside_canonical_transcript", details

    overlapping_exon_indexes = []
    for exon_index, exon in enumerate(exons):
        if exon["start_1based"] <= dipa_position <= exon["end_1based"]:
            overlapping_exon_indexes.append(exon_index)

    removed_exon = None

    if len(overlapping_exon_indexes) > 1:
        return None, "dipa_overlaps_multiple_canonical_exons", ""

    if len(overlapping_exon_indexes) == 1:
        removed_exon_index = overlapping_exon_indexes[0]
        removed_exon = exons[removed_exon_index]
        upstream_exons = exons[:removed_exon_index]
        downstream_exons = exons[removed_exon_index + 1 :]
        dipa_context = "exonic_exon_removed"
    else:
        upstream_exons = None
        downstream_exons = None

        for exon_index in range(len(exons) - 1):
            first_exon = exons[exon_index]
            second_exon = exons[exon_index + 1]

            # This genomic gap works for both strands. The exon list itself is
            # already in transcript order.
            intron_start = min(
                first_exon["end_1based"],
                second_exon["end_1based"],
            ) + 1
            intron_end = max(
                first_exon["start_1based"],
                second_exon["start_1based"],
            ) - 1

            if intron_start <= dipa_position <= intron_end:
                upstream_exons = exons[: exon_index + 1]
                downstream_exons = exons[exon_index + 1 :]
                break

        if upstream_exons is None or downstream_exons is None:
            return None, "dipa_not_in_canonical_exon_or_intron", ""

        dipa_context = "intronic"

    upstream_exonic_bp = sum(
        exon["end_1based"] - exon["start_1based"] + 1
        for exon in upstream_exons
    )
    downstream_exonic_bp = sum(
        exon["end_1based"] - exon["start_1based"] + 1
        for exon in downstream_exons
    )

    upstream_is_short = upstream_exonic_bp < min_side_bp
    downstream_is_short = downstream_exonic_bp < min_side_bp

    if upstream_is_short and downstream_is_short:
        details = (
            f"upstream={upstream_exonic_bp} bp; "
            f"downstream={downstream_exonic_bp} bp"
        )
        return None, "insufficient_cdna_on_both_sides", details

    if upstream_is_short:
        details = f"upstream={upstream_exonic_bp} bp"
        return None, "insufficient_upstream_cdna", details

    if downstream_is_short:
        details = f"downstream={downstream_exonic_bp} bp"
        return None, "insufficient_downstream_cdna", details

    transcript_gene_ids = transcript.attributes.get("gene_id", [])
    transcript_gene_id = (
        transcript_gene_ids[0] if transcript_gene_ids else ""
    )

    removed_exon_id = ""
    removed_exon_start = ""
    removed_exon_end = ""
    if removed_exon is not None:
        removed_exon_id = removed_exon["exon_id"]
        removed_exon_start = removed_exon["start_1based"]
        removed_exon_end = removed_exon["end_1based"]

    legacy_pas_strand = apa_row["_legacy_pas_strand"]
    if pd.isna(legacy_pas_strand) or legacy_pas_strand not in {"+", "-"}:
        legacy_pas_strand = ""
        legacy_strand_matches_annotation = "unavailable"
    elif legacy_pas_strand == transcript.strand:
        legacy_strand_matches_annotation = "yes"
    else:
        legacy_strand_matches_annotation = "no"

    gene_model = {
        "gene_symbol": str(apa_row["gene_symbol"]).strip(),
        "gene_id": transcript_gene_id,
        "transcript_id": transcript.id,
        "chromosome": transcript.seqid,
        "strand": transcript.strand,
        "strand_source": "Ensembl_canonical_transcript",
        "legacy_pas_strand": legacy_pas_strand,
        "legacy_strand_matches_annotation": (
            legacy_strand_matches_annotation
        ),
        "PASid": apa_row["PASid"],
        "RED": float(apa_row["_RED_numeric"]),
        "p_adj": float(apa_row["_p_adj_numeric"]),
        "dipa_start_0based": int(apa_row["_start_numeric"]),
        "dipa_position_1based": dipa_position,
        "dipa_context": dipa_context,
        "removed_exon_id": removed_exon_id,
        "removed_exon_start": removed_exon_start,
        "removed_exon_end": removed_exon_end,
        "upstream_exonic_bp": upstream_exonic_bp,
        "downstream_exonic_bp": downstream_exonic_bp,
        "upstream_exons": upstream_exons,
        "downstream_exons": downstream_exons,
        "source_row": {
            key: value
            for key, value in apa_row.to_dict().items()
            if not str(key).startswith("_")
        },
    }

    return gene_model, None, ""


def inspect_bigwig(sample_record):
    """Open one BigWig and return chromosome lengths for preflight checks."""

    bigwig_path = sample_record["_bigwig_path"]
    try:
        with pyBigWig.open(bigwig_path) as bigwig_file:
            chromosome_lengths = bigwig_file.chroms()
    except Exception as error:
        raise RuntimeError(
            f"Could not open BigWig for sample "
            f"{sample_record['sample_id']}: {bigwig_path}. {error}"
        ) from error

    if not chromosome_lengths:
        raise RuntimeError(
            f"BigWig contains no chromosome information: {bigwig_path}"
        )

    return sample_record["sample_id"], chromosome_lengths


def bin_cdna_signal(cdna_values, bins_per_side):
    """Average a cDNA signal array into equal or nearly equal nonempty bins."""

    if len(cdna_values) < bins_per_side:
        raise ValueError(
            f"Cannot divide {len(cdna_values)} bases into "
            f"{bins_per_side} nonempty bins."
        )

    signal_groups = np.array_split(cdna_values, bins_per_side)
    return np.array(
        [float(np.mean(signal_group)) for signal_group in signal_groups],
        dtype=float,
    )


def extract_oriented_exon_values(
    bigwig_file,
    sample_id,
    gene_model,
    exons,
):
    """Extract, validate, and orient transcript exon signal."""

    signal_parts = []

    for exon in exons:
        exon_values = bigwig_file.values(
            gene_model["chromosome"],
            exon["start_0based"],
            exon["end_0based"],
            numpy=True,
        )
        exon_values = np.asarray(exon_values, dtype=float)
        exon_values = np.nan_to_num(exon_values, nan=0.0)

        if np.any(exon_values < 0):
            raise ValueError(
                f"Negative BigWig value found for sample {sample_id}, "
                f"gene {gene_model['gene_symbol']}, exon "
                f"{exon['exon_id']}. Log2-ratio normalization requires "
                f"nonnegative linear signal."
            )

        if gene_model["strand"] == "-":
            exon_values = exon_values[::-1]

        signal_parts.append(exon_values)

    return np.concatenate(signal_parts)


def extract_sample_profiles(sample_record, gene_models, bins_per_side):
    """
    Extract and bin every retained gene for one BigWig sample.

    Each thread calls this function for one sample and opens its own BigWig
    handle. BigWig handles are never shared between threads.
    """

    sample_id = sample_record["sample_id"]
    bigwig_path = sample_record["_bigwig_path"]
    total_bins = bins_per_side * 2
    sample_profiles = np.empty(
        (len(gene_models), total_bins),
        dtype=float,
    )

    try:
        bigwig_file = pyBigWig.open(bigwig_path)
    except Exception as error:
        raise RuntimeError(
            f"Could not open BigWig for sample {sample_id}: "
            f"{bigwig_path}. {error}"
        ) from error

    if bigwig_file is None:
        raise RuntimeError(
            f"pyBigWig could not open the file for sample {sample_id}: "
            f"{bigwig_path}"
        )

    try:
        for gene_index, gene_model in enumerate(gene_models):
            upstream_signal = extract_oriented_exon_values(
                bigwig_file,
                sample_id,
                gene_model,
                gene_model["upstream_exons"],
            )
            downstream_signal = extract_oriented_exon_values(
                bigwig_file,
                sample_id,
                gene_model,
                gene_model["downstream_exons"],
            )

            upstream_bins = bin_cdna_signal(
                upstream_signal,
                bins_per_side,
            )
            downstream_bins = bin_cdna_signal(
                downstream_signal,
                bins_per_side,
            )

            sample_profiles[gene_index, :] = np.concatenate(
                [upstream_bins, downstream_bins]
            )
    finally:
        bigwig_file.close()

    return sample_id, sample_profiles


def choose_thread_count(requested_threads, number_of_samples):
    """Choose an HPC-friendly thread count without exceeding sample count."""

    if requested_threads < 0:
        raise ValueError("--threads must be zero or a positive integer.")

    if requested_threads > 0:
        return min(requested_threads, number_of_samples)

    scheduler_variables = [
        "SLURM_CPUS_PER_TASK",
        "PBS_NP",
        "NSLOTS",
    ]

    detected_cpus = None
    for variable_name in scheduler_variables:
        variable_value = os.environ.get(variable_name)
        if variable_value:
            try:
                detected_cpus = int(variable_value)
            except ValueError:
                continue
            if detected_cpus > 0:
                break

    if detected_cpus is None:
        detected_cpus = os.cpu_count() or 1

    return max(1, min(detected_cpus, number_of_samples))


def write_per_gene_log2fc(
    output_path,
    gene_models,
    log2fc_matrix,
    treatment_record,
    position_indexes,
    regions,
    normalized_positions,
    write_header,
):
    """Write gene-by-position output in memory-bounded chunks."""

    genes_per_output_chunk = 500
    total_positions = len(position_indexes)

    for first_gene in range(0, len(gene_models), genes_per_output_chunk):
        last_gene = min(
            first_gene + genes_per_output_chunk,
            len(gene_models),
        )
        chunk_models = gene_models[first_gene:last_gene]
        chunk_matrix = log2fc_matrix[first_gene:last_gene, :]
        chunk_gene_symbols = [
            gene_model["gene_symbol"] for gene_model in chunk_models
        ]

        output_frame = pd.DataFrame(
            {
                "gene_symbol": np.repeat(
                    chunk_gene_symbols,
                    total_positions,
                ),
                "condition": treatment_record["condition"],
                "sample_id": treatment_record["sample_id"],
                "pair_id": treatment_record["pair_id"],
                "position_index": np.tile(
                    position_indexes,
                    len(chunk_models),
                ),
                "region": np.tile(regions, len(chunk_models)),
                "normalized_position": np.tile(
                    normalized_positions,
                    len(chunk_models),
                ),
                "log2fc": chunk_matrix.reshape(-1),
            }
        )

        output_frame.to_csv(
            output_path,
            sep="\t",
            index=False,
            mode="w" if write_header else "a",
            header=write_header,
        )
        write_header = False

    return write_header


def get_package_versions():
    """Collect versions for reproducibility in run_parameters.json."""

    package_names = [
        "numpy",
        "pandas",
        "gffutils",
        "pyBigWig",
        "matplotlib",
    ]
    versions = {}

    for package_name in package_names:
        try:
            versions[package_name] = version(package_name)
        except PackageNotFoundError:
            versions[package_name] = "unknown"

    versions["python"] = sys.version.split()[0]
    return versions


def parse_arguments():
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Create exon-only metagene profiles around significant UP dIPAs "
            "using per-replicate BigWig signal."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--apa-results",
        required=True,
        help="Tab-separated APAlyzer results table.",
    )
    parser.add_argument(
        "--gtf",
        required=True,
        help="GENCODE GTF using the same genome assembly as the BigWigs.",
    )
    parser.add_argument(
        "--samples",
        required=True,
        help="Tab-separated sample sheet.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for plots, tables, QC, and the reusable GTF database.",
    )
    parser.add_argument(
        "--pseudocount",
        required=True,
        type=float,
        help="Positive value added before calculating BigWig log2 ratios.",
    )
    parser.add_argument(
        "--adjusted-pvalue",
        type=float,
        default=0.05,
        help="Maximum APAlyzer adjusted p-value.",
    )
    parser.add_argument(
        "--bins-per-side",
        type=int,
        default=100,
        help="Number of normalized bins before and after the dIPA.",
    )
    parser.add_argument(
        "--min-side-bp",
        type=int,
        default=100,
        help="Minimum retained exonic bases on each side of the dIPA.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=0,
        help=(
            "Number of sample-level extraction threads. Zero detects common "
            "HPC scheduler CPU variables and otherwise uses available CPUs."
        ),
    )
    parser.add_argument(
        "--gtf-db",
        default="",
        help=(
            "Optional path for the reusable gffutils SQLite database. The "
            "default places it inside the output directory."
        ),
    )
    parser.add_argument(
        "--rebuild-gtf-db",
        action="store_true",
        help="Delete and rebuild an existing gffutils database.",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    if args.pseudocount <= 0:
        raise ValueError("--pseudocount must be greater than zero.")
    if not 0 < args.adjusted_pvalue <= 1:
        raise ValueError("--adjusted-pvalue must be greater than 0 and at most 1.")
    if args.bins_per_side <= 0:
        raise ValueError("--bins-per-side must be a positive integer.")
    if args.min_side_bp <= 0:
        raise ValueError("--min-side-bp must be a positive integer.")
    if args.min_side_bp < args.bins_per_side:
        raise ValueError(
            "--min-side-bp must be at least --bins-per-side so every "
            "normalized bin contains at least one real nucleotide."
        )

    apa_results_path = Path(args.apa_results).expanduser().resolve()
    gtf_path = Path(args.gtf).expanduser().resolve()
    sample_sheet_path = Path(args.samples).expanduser().resolve()
    output_directory = Path(args.output_dir).expanduser().resolve()

    for input_path, input_label in [
        (apa_results_path, "APAlyzer results"),
        (gtf_path, "GTF"),
        (sample_sheet_path, "sample sheet"),
    ]:
        if not input_path.is_file():
            raise FileNotFoundError(
                f"{input_label} file does not exist: {input_path}"
            )

    output_directory.mkdir(parents=True, exist_ok=True)

    if args.gtf_db:
        gtf_database_path = Path(args.gtf_db).expanduser().resolve()
        gtf_database_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        gtf_database_path = output_directory / (
            gtf_path.name + ".gffutils.db"
        )

    print("Reading APAlyzer results...")
    apa_results = pd.read_csv(
        apa_results_path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    missing_apa_columns = [
        column
        for column in REQUIRED_APA_COLUMNS
        if column not in apa_results.columns
    ]
    if missing_apa_columns:
        raise ValueError(
            "APAlyzer results are missing required columns: "
            + ", ".join(missing_apa_columns)
        )

    numeric_apa_columns = {
        "RED": "_RED_numeric",
        "p_adj": "_p_adj_numeric",
        "start": "_start_numeric",
        "end": "_end_numeric",
        "browser_start_1based": "_browser_start_numeric",
    }

    for input_column, helper_column in numeric_apa_columns.items():
        try:
            apa_results[helper_column] = pd.to_numeric(
                apa_results[input_column],
                errors="raise",
            )
        except Exception as error:
            raise ValueError(
                f"APAlyzer column {input_column} contains a nonnumeric value."
            ) from error

    for integer_helper_column in [
        "_start_numeric",
        "_end_numeric",
        "_browser_start_numeric",
    ]:
        integer_values = apa_results[integer_helper_column]
        if not np.all(np.equal(integer_values, np.floor(integer_values))):
            raise ValueError(
                f"APAlyzer coordinate column {integer_helper_column} "
                "contains a non-integer value."
            )
        apa_results[integer_helper_column] = integer_values.astype(np.int64)

    apa_results["_APAreg_clean"] = (
        apa_results["APAreg"].astype(str).str.strip().str.upper()
    )
    apa_results["_APAreg_original_clean"] = (
        apa_results["APAreg_original"]
        .astype(str)
        .str.strip()
        .str.upper()
    )
    apa_results["_pas_match_boolean"] = (
        apa_results["pas_coordinate_match"]
        .astype(str)
        .str.strip()
        .str.lower()
        .isin(TRUE_TEXT_VALUES)
    )
    # PASid may be a label from an older genome assembly. Its strand is kept
    # only for QC and is never used to select or orient the current transcript.
    apa_results["_legacy_pas_strand"] = apa_results["PASid"].apply(
        parse_pas_strand
    )

    excluded_rows = []
    qualifying_indexes = []

    for row_index, row in apa_results.iterrows():
        if float(row["_p_adj_numeric"]) > args.adjusted_pvalue:
            add_exclusion(
                excluded_rows,
                row,
                "APA_filtering",
                "not_significant",
                f"p_adj={row['_p_adj_numeric']}",
            )
            continue

        if row["_APAreg_clean"] != "UP":
            add_exclusion(
                excluded_rows,
                row,
                "APA_filtering",
                "not_UP",
                f"APAreg={row['APAreg']}",
            )
            continue

        if float(row["_RED_numeric"]) <= 0:
            add_exclusion(
                excluded_rows,
                row,
                "APA_filtering",
                "nonpositive_RED_for_UP",
                f"RED={row['_RED_numeric']}",
            )
            continue

        if not bool(row["_pas_match_boolean"]):
            add_exclusion(
                excluded_rows,
                row,
                "APA_filtering",
                "pas_coordinate_match_false",
                f"pas_coordinate_match={row['pas_coordinate_match']}",
            )
            continue

        if row["_APAreg_original_clean"] != "UP":
            add_exclusion(
                excluded_rows,
                row,
                "APA_filtering",
                "APAreg_and_APAreg_original_disagree",
                (
                    f"APAreg={row['APAreg']}; "
                    f"APAreg_original={row['APAreg_original']}"
                ),
            )
            continue

        if not str(row["gene_symbol"]).strip():
            add_exclusion(
                excluded_rows,
                row,
                "APA_filtering",
                "missing_gene_symbol",
            )
            continue

        if int(row["_end_numeric"]) != int(row["_start_numeric"]) + 1:
            add_exclusion(
                excluded_rows,
                row,
                "coordinate_validation",
                "dipa_interval_is_not_one_base",
                f"start={row['start']}; end={row['end']}",
            )
            continue

        if (
            int(row["_browser_start_numeric"])
            != int(row["_start_numeric"]) + 1
        ):
            add_exclusion(
                excluded_rows,
                row,
                "coordinate_validation",
                "bed_and_gtf_coordinates_disagree",
                (
                    f"start={row['start']}; "
                    f"browser_start_1based={row['browser_start_1based']}"
                ),
            )
            continue

        qualifying_indexes.append(row_index)

    if not qualifying_indexes:
        excluded_output_path = output_directory / "excluded_genes.tsv"
        pd.DataFrame(excluded_rows).to_csv(
            excluded_output_path,
            sep="\t",
            index=False,
        )
        raise ValueError(
            "No significant UP dIPAs remained after APA and coordinate "
            f"filtering. QC was written to {excluded_output_path}"
        )

    qualifying_rows = apa_results.loc[qualifying_indexes].copy()
    qualifying_rows["_gene_symbol_clean"] = (
        qualifying_rows["gene_symbol"].astype(str).str.strip()
    )
    qualifying_rows = qualifying_rows.sort_values(
        by=[
            "_gene_symbol_clean",
            "_RED_numeric",
            "_p_adj_numeric",
            "_start_numeric",
        ],
        ascending=[True, False, True, True],
        kind="mergesort",
    )

    selected_apa_indexes = []
    for gene_symbol, gene_rows in qualifying_rows.groupby(
        "_gene_symbol_clean",
        sort=False,
    ):
        selected_apa_indexes.append(gene_rows.index[0])

        for additional_index in gene_rows.index[1:]:
            add_exclusion(
                excluded_rows,
                qualifying_rows.loc[additional_index],
                "representative_dipa_selection",
                "additional_significant_PAS_for_gene",
                (
                    f"gene={gene_symbol}; selected highest RED site "
                    f"{gene_rows.iloc[0]['PASid']}"
                ),
            )

    selected_apa_rows = qualifying_rows.loc[selected_apa_indexes].copy()

    gtf_database = create_or_open_gtf_database(
        gtf_path,
        gtf_database_path,
        args.rebuild_gtf_db,
    )

    print("Indexing Ensembl_canonical transcripts...")
    canonical_by_gene_id = {}
    canonical_by_symbol_and_chromosome = {}

    for transcript in gtf_database.features_of_type("transcript"):
        transcript_tags = transcript.attributes.get("tag", [])
        if "Ensembl_canonical" not in transcript_tags:
            continue

        transcript_gene_ids = transcript.attributes.get("gene_id", [])
        transcript_gene_names = transcript.attributes.get("gene_name", [])

        if transcript_gene_ids:
            normalized_id = normalize_gene_id(transcript_gene_ids[0])
            canonical_by_gene_id.setdefault(normalized_id, []).append(
                transcript
            )

        if transcript_gene_names:
            symbol_key = (
                transcript_gene_names[0],
                transcript.seqid,
            )
            canonical_by_symbol_and_chromosome.setdefault(
                symbol_key,
                [],
            ).append(transcript)

    gene_models = []

    for _, apa_row in selected_apa_rows.iterrows():
        candidate_transcripts = []
        input_gene_id = ""

        if "gene_id" in selected_apa_rows.columns:
            input_gene_id = str(apa_row["gene_id"]).strip()
            if input_gene_id.lower() in {"", ".", "na", "nan", "none"}:
                input_gene_id = ""

        if input_gene_id:
            candidate_transcripts = canonical_by_gene_id.get(
                normalize_gene_id(input_gene_id),
                [],
            )
            candidate_transcripts = [
                transcript
                for transcript in candidate_transcripts
                if transcript.seqid == str(apa_row["chr"]).strip()
            ]
        else:
            symbol_key = (
                str(apa_row["gene_symbol"]).strip(),
                str(apa_row["chr"]).strip(),
            )
            candidate_transcripts = (
                canonical_by_symbol_and_chromosome.get(symbol_key, [])
            )

        # The mapped coordinate and current annotation are authoritative. This
        # span check also disambiguates same-name genes on opposite strands.
        dipa_position = int(apa_row["_browser_start_numeric"])
        candidate_transcripts = [
            transcript
            for transcript in candidate_transcripts
            if int(transcript.start) <= dipa_position <= int(transcript.end)
        ]

        if not candidate_transcripts:
            add_exclusion(
                excluded_rows,
                apa_row,
                "canonical_transcript_selection",
                "no_matching_Ensembl_canonical_transcript",
                (
                    f"gene={apa_row['gene_symbol']}; chr={apa_row['chr']}; "
                    f"mapped_position={dipa_position}"
                ),
            )
            continue

        if len(candidate_transcripts) > 1:
            transcript_ids = ", ".join(
                (
                    f"{transcript.id}({transcript.strand}:"
                    f"{transcript.start}-{transcript.end})"
                )
                for transcript in candidate_transcripts
            )
            add_exclusion(
                excluded_rows,
                apa_row,
                "canonical_transcript_selection",
                "multiple_matching_Ensembl_canonical_transcripts",
                transcript_ids,
            )
            continue

        canonical_transcript = candidate_transcripts[0]
        gene_model, exclusion_reason, exclusion_details = prepare_cdna_model(
            apa_row,
            canonical_transcript,
            gtf_database,
            args.min_side_bp,
        )

        if exclusion_reason is not None:
            add_exclusion(
                excluded_rows,
                apa_row,
                "cdna_model_construction",
                exclusion_reason,
                exclusion_details,
            )
            continue

        gene_models.append(gene_model)

    if not gene_models:
        excluded_output_path = output_directory / "excluded_genes.tsv"
        pd.DataFrame(excluded_rows).to_csv(
            excluded_output_path,
            sep="\t",
            index=False,
        )
        raise ValueError(
            "No genes remained after canonical transcript and cDNA filtering. "
            f"QC was written to {excluded_output_path}"
        )

    print("Reading sample sheet...")
    sample_sheet = pd.read_csv(
        sample_sheet_path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    missing_sample_columns = [
        column
        for column in REQUIRED_SAMPLE_COLUMNS
        if column not in sample_sheet.columns
    ]
    if missing_sample_columns:
        raise ValueError(
            "Sample sheet is missing required columns: "
            + ", ".join(missing_sample_columns)
        )

    for text_column in REQUIRED_SAMPLE_COLUMNS:
        sample_sheet[text_column] = (
            sample_sheet[text_column].astype(str).str.strip()
        )

    if sample_sheet["sample_id"].eq("").any():
        raise ValueError("Sample sheet contains a blank sample_id.")
    if sample_sheet["sample_id"].duplicated().any():
        duplicated_ids = sample_sheet.loc[
            sample_sheet["sample_id"].duplicated(keep=False),
            "sample_id",
        ].unique()
        raise ValueError(
            "Sample sheet contains duplicate sample_id values: "
            + ", ".join(duplicated_ids)
        )
    if sample_sheet["condition"].eq("").any():
        raise ValueError("Sample sheet contains a blank condition.")

    sample_sheet["role"] = sample_sheet["role"].str.lower()
    invalid_roles = sorted(
        set(sample_sheet["role"]) - {"control", "treatment"}
    )
    if invalid_roles:
        raise ValueError(
            "Sample sheet role must be control or treatment. Invalid values: "
            + ", ".join(invalid_roles)
        )

    control_samples = sample_sheet[sample_sheet["role"] == "control"]
    treatment_samples = sample_sheet[sample_sheet["role"] == "treatment"]

    if control_samples.empty:
        raise ValueError("Sample sheet contains no control samples.")
    if treatment_samples.empty:
        raise ValueError("Sample sheet contains no treatment samples.")

    sample_sheet_directory = sample_sheet_path.parent
    resolved_bigwig_paths = []
    for bigwig_text in sample_sheet["bigwig"]:
        bigwig_path = Path(bigwig_text).expanduser()
        if not bigwig_path.is_absolute():
            bigwig_path = sample_sheet_directory / bigwig_path
        bigwig_path = bigwig_path.resolve()
        if not bigwig_path.is_file():
            raise FileNotFoundError(
                f"BigWig file does not exist: {bigwig_path}"
            )
        resolved_bigwig_paths.append(str(bigwig_path))

    sample_sheet["_bigwig_path"] = resolved_bigwig_paths
    control_samples = sample_sheet[sample_sheet["role"] == "control"]
    treatment_samples = sample_sheet[sample_sheet["role"] == "treatment"]

    treatment_has_pair = treatment_samples["pair_id"].ne("")
    if treatment_has_pair.any() and not treatment_has_pair.all():
        raise ValueError(
            "Do not mix paired and unpaired treatment samples. Fill pair_id "
            "for every treatment or leave every pair_id blank."
        )

    paired_controls = bool(treatment_has_pair.all())

    if paired_controls:
        if control_samples["pair_id"].eq("").any():
            raise ValueError(
                "Paired mode requires pair_id for every control sample."
            )

        control_pair_counts = control_samples.groupby("pair_id").size()
        duplicated_control_pairs = control_pair_counts[
            control_pair_counts != 1
        ]
        if not duplicated_control_pairs.empty:
            raise ValueError(
                "Paired mode requires exactly one control for each pair_id. "
                "Problem pair_id values: "
                + ", ".join(duplicated_control_pairs.index.astype(str))
            )

        available_control_pairs = set(control_samples["pair_id"])
        missing_control_pairs = sorted(
            set(treatment_samples["pair_id"]) - available_control_pairs
        )
        if missing_control_pairs:
            raise ValueError(
                "No matching control exists for treatment pair_id values: "
                + ", ".join(missing_control_pairs)
            )
    else:
        if sample_sheet["pair_id"].ne("").any():
            raise ValueError(
                "Unpaired mode requires pair_id to be blank for every sample."
            )
        print(
            "WARNING: Controls are unpaired. Each treatment replicate will be "
            "compared with the mean control profile, so the plotted SD does "
            "not fully preserve control-replicate variation."
        )

    sample_records = sample_sheet.to_dict(orient="records")
    thread_count = choose_thread_count(args.threads, len(sample_records))
    print(
        f"Using {thread_count} extraction thread(s) for "
        f"{len(sample_records)} BigWig sample(s)."
    )

    print("Checking BigWig chromosome compatibility...")
    bigwig_chromosomes = {}
    with ThreadPoolExecutor(max_workers=thread_count) as thread_pool:
        inspection_futures = {
            thread_pool.submit(inspect_bigwig, sample_record): sample_record
            for sample_record in sample_records
        }
        for completed_future in as_completed(inspection_futures):
            sample_id, chromosome_lengths = completed_future.result()
            bigwig_chromosomes[sample_id] = chromosome_lengths

    compatible_gene_models = []
    for gene_model in gene_models:
        compatibility_problem = None

        all_gene_exons = (
            gene_model["upstream_exons"]
            + gene_model["downstream_exons"]
        )
        maximum_exon_end = max(
            exon["end_0based"] for exon in all_gene_exons
        )

        for sample_record in sample_records:
            sample_id = sample_record["sample_id"]
            chromosome_lengths = bigwig_chromosomes[sample_id]
            chromosome = gene_model["chromosome"]

            if chromosome not in chromosome_lengths:
                compatibility_problem = (
                    "bigwig_chromosome_missing",
                    f"sample={sample_id}; chromosome={chromosome}",
                )
                break

            if maximum_exon_end > chromosome_lengths[chromosome]:
                compatibility_problem = (
                    "bigwig_coordinate_out_of_bounds",
                    (
                        f"sample={sample_id}; chromosome={chromosome}; "
                        f"required_end={maximum_exon_end}; "
                        f"chromosome_length={chromosome_lengths[chromosome]}"
                    ),
                )
                break

        if compatibility_problem is not None:
            add_exclusion(
                excluded_rows,
                gene_model["source_row"],
                "bigwig_preflight",
                compatibility_problem[0],
                compatibility_problem[1],
            )
            continue

        compatible_gene_models.append(gene_model)

    gene_models = compatible_gene_models

    if not gene_models:
        excluded_output_path = output_directory / "excluded_genes.tsv"
        pd.DataFrame(excluded_rows).to_csv(
            excluded_output_path,
            sep="\t",
            index=False,
        )
        raise ValueError(
            "No genes remained after BigWig chromosome checks. "
            f"QC was written to {excluded_output_path}"
        )

    print(
        f"Extracting BigWig signal for {len(gene_models)} genes "
        f"and {len(sample_records)} samples..."
    )
    sample_signal_matrices = {}

    with ThreadPoolExecutor(max_workers=thread_count) as thread_pool:
        extraction_futures = {
            thread_pool.submit(
                extract_sample_profiles,
                sample_record,
                gene_models,
                args.bins_per_side,
            ): sample_record
            for sample_record in sample_records
        }

        for completed_future in as_completed(extraction_futures):
            sample_record = extraction_futures[completed_future]
            try:
                sample_id, sample_profiles = completed_future.result()
            except Exception as error:
                raise RuntimeError(
                    f"Signal extraction failed for sample "
                    f"{sample_record['sample_id']}: {error}"
                ) from error

            sample_signal_matrices[sample_id] = sample_profiles
            print(f"  Finished {sample_id}")

    selected_output_columns = [
        "gene_symbol",
        "gene_id",
        "transcript_id",
        "chromosome",
        "strand",
        "strand_source",
        "legacy_pas_strand",
        "legacy_strand_matches_annotation",
        "PASid",
        "RED",
        "p_adj",
        "dipa_start_0based",
        "dipa_position_1based",
        "dipa_context",
        "removed_exon_id",
        "removed_exon_start",
        "removed_exon_end",
        "upstream_exonic_bp",
        "downstream_exonic_bp",
    ]
    selected_output_rows = [
        {
            output_column: gene_model[output_column]
            for output_column in selected_output_columns
        }
        for gene_model in gene_models
    ]
    pd.DataFrame(selected_output_rows).to_csv(
        output_directory / "selected_up_genes.tsv",
        sep="\t",
        index=False,
    )

    if excluded_rows:
        excluded_frame = pd.DataFrame(excluded_rows)
    else:
        excluded_frame = pd.DataFrame(
            columns=REQUIRED_APA_COLUMNS
            + ["exclusion_stage", "exclusion_reason", "details"]
        )
    excluded_frame.to_csv(
        output_directory / "excluded_genes.tsv",
        sep="\t",
        index=False,
    )

    total_positions = args.bins_per_side * 2
    position_indexes = np.arange(total_positions)
    regions = np.where(
        position_indexes < args.bins_per_side,
        "upstream",
        "downstream",
    )
    normalized_positions = (
        position_indexes.astype(float) + 0.5
    ) / total_positions

    control_by_pair = {}
    if paired_controls:
        for control_record in control_samples.to_dict(orient="records"):
            control_by_pair[control_record["pair_id"]] = control_record[
                "sample_id"
            ]
        mean_control_matrix = None
    else:
        control_matrices = [
            sample_signal_matrices[sample_id]
            for sample_id in control_samples["sample_id"]
        ]
        mean_control_matrix = np.mean(
            np.stack(control_matrices, axis=0),
            axis=0,
        )

    per_gene_output_path = output_directory / "per_gene_log2fc.tsv"
    write_per_gene_header = True
    replicate_output_frames = []
    condition_replicate_profiles = {}

    treatment_records = treatment_samples.to_dict(orient="records")
    for treatment_record in treatment_records:
        treatment_matrix = sample_signal_matrices[
            treatment_record["sample_id"]
        ]

        if paired_controls:
            control_sample_id = control_by_pair[treatment_record["pair_id"]]
            control_matrix = sample_signal_matrices[control_sample_id]
        else:
            control_matrix = mean_control_matrix

        log2fc_matrix = np.log2(
            (treatment_matrix + args.pseudocount)
            / (control_matrix + args.pseudocount)
        )

        write_per_gene_header = write_per_gene_log2fc(
            per_gene_output_path,
            gene_models,
            log2fc_matrix,
            treatment_record,
            position_indexes,
            regions,
            normalized_positions,
            write_per_gene_header,
        )

        replicate_profile = np.mean(log2fc_matrix, axis=0)
        condition_name = treatment_record["condition"]
        condition_replicate_profiles.setdefault(
            condition_name,
            [],
        ).append(replicate_profile)

        replicate_output_frames.append(
            pd.DataFrame(
                {
                    "condition": condition_name,
                    "sample_id": treatment_record["sample_id"],
                    "pair_id": treatment_record["pair_id"],
                    "position_index": position_indexes,
                    "region": regions,
                    "normalized_position": normalized_positions,
                    "mean_gene_log2fc": replicate_profile,
                    "number_of_genes": len(gene_models),
                }
            )
        )

    replicate_profiles = pd.concat(
        replicate_output_frames,
        ignore_index=True,
    )
    replicate_profiles.to_csv(
        output_directory / "replicate_profiles.tsv",
        sep="\t",
        index=False,
    )

    treatment_summary_frames = []
    treatment_condition_order = list(
        treatment_samples["condition"].drop_duplicates()
    )

    for condition_name in treatment_condition_order:
        condition_matrix = np.stack(
            condition_replicate_profiles[condition_name],
            axis=0,
        )
        mean_profile = np.mean(condition_matrix, axis=0)

        if condition_matrix.shape[0] > 1:
            standard_deviation = np.std(
                condition_matrix,
                axis=0,
                ddof=1,
            )
        else:
            standard_deviation = np.full(total_positions, np.nan)

        treatment_summary_frames.append(
            pd.DataFrame(
                {
                    "condition": condition_name,
                    "position_index": position_indexes,
                    "region": regions,
                    "normalized_position": normalized_positions,
                    "mean_log2fc": mean_profile,
                    "sd_log2fc": standard_deviation,
                    "number_of_replicates": condition_matrix.shape[0],
                    "number_of_genes": len(gene_models),
                }
            )
        )

    treatment_summary = pd.concat(
        treatment_summary_frames,
        ignore_index=True,
    )
    treatment_summary.to_csv(
        output_directory / "treatment_summary.tsv",
        sep="\t",
        index=False,
    )

    print("Drawing combined metagene plot...")
    figure, axis = plt.subplots(figsize=(9, 5.5))

    for condition_index, condition_name in enumerate(
        treatment_condition_order
    ):
        condition_summary = treatment_summary[
            treatment_summary["condition"] == condition_name
        ]
        color = COLORBLIND_COLORS[
            condition_index % len(COLORBLIND_COLORS)
        ]
        x_values = condition_summary["normalized_position"].to_numpy()
        mean_values = condition_summary["mean_log2fc"].to_numpy()
        sd_values = condition_summary["sd_log2fc"].to_numpy()

        axis.plot(
            x_values,
            mean_values,
            color=color,
            linewidth=2,
            label=condition_name,
        )

        if not np.all(np.isnan(sd_values)):
            axis.fill_between(
                x_values,
                mean_values - sd_values,
                mean_values + sd_values,
                color=color,
                alpha=0.22,
                linewidth=0,
            )

    axis.axhline(0, color="#555555", linewidth=1, linestyle="--")
    axis.axvline(0.5, color="#333333", linewidth=1, linestyle=":")
    axis.set_xlim(0, 1)
    axis.set_xticks([0, 0.5, 1])
    axis.set_xticklabels(["0%", "dIPA", "100%"])
    axis.set_xlabel("Normalized canonical cDNA position")
    axis.set_ylabel("log2 fold-change vs control")
    axis.set_title("Metagene signal around UP-regulated dIPAs")
    axis.legend(frameon=False)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.grid(axis="y", color="#DDDDDD", linewidth=0.6)

    replicate_counts = treatment_samples.groupby("condition").size()
    replicate_text = ", ".join(
        f"{condition}: {replicate_counts[condition]} replicate(s)"
        for condition in treatment_condition_order
    )
    figure.text(
        0.5,
        0.01,
        f"{len(gene_models)} genes; {replicate_text}; shading = +/- 1 SD",
        ha="center",
        fontsize=9,
    )
    figure.tight_layout(rect=[0, 0.04, 1, 1])
    figure.savefig(
        output_directory / "combined_metagene.png",
        dpi=300,
    )
    figure.savefig(output_directory / "combined_metagene.pdf")
    plt.close(figure)

    intronic_count = sum(
        gene_model["dipa_context"] == "intronic"
        for gene_model in gene_models
    )
    exonic_count = sum(
        gene_model["dipa_context"] == "exonic_exon_removed"
        for gene_model in gene_models
    )
    legacy_strand_mismatch_count = sum(
        gene_model["legacy_strand_matches_annotation"] == "no"
        for gene_model in gene_models
    )
    legacy_strand_unavailable_count = sum(
        gene_model["legacy_strand_matches_annotation"] == "unavailable"
        for gene_model in gene_models
    )

    filtering_counts = {
        "apa_rows_read": int(len(apa_results)),
        "significant_up_pas_rows": int(len(qualifying_rows)),
        "unique_significant_up_genes": int(
            qualifying_rows["_gene_symbol_clean"].nunique()
        ),
        "final_genes": int(len(gene_models)),
        "final_intronic_dipas": int(intronic_count),
        "final_exonic_dipas": int(exonic_count),
        "legacy_pas_strand_mismatches": int(
            legacy_strand_mismatch_count
        ),
        "legacy_pas_strand_unavailable": int(
            legacy_strand_unavailable_count
        ),
    }

    run_parameters = {
        "run_time_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "apa_results": str(apa_results_path),
            "gtf": str(gtf_path),
            "samples": str(sample_sheet_path),
            "gtf_database": str(gtf_database_path),
        },
        "output_directory": str(output_directory),
        "parameters": {
            "pseudocount": args.pseudocount,
            "adjusted_pvalue": args.adjusted_pvalue,
            "bins_per_side": args.bins_per_side,
            "min_side_bp": args.min_side_bp,
            "requested_threads": args.threads,
            "threads_used": thread_count,
            "paired_controls": paired_controls,
        },
        "filtering_counts": filtering_counts,
        "treatments": {
            condition: int(replicate_counts[condition])
            for condition in treatment_condition_order
        },
        "package_versions": get_package_versions(),
    }

    with (output_directory / "run_parameters.json").open(
        "w",
        encoding="utf-8",
    ) as output_file:
        json.dump(run_parameters, output_file, indent=2)

    short_cdna_exclusions = sum(
        row["exclusion_reason"]
        in {
            "insufficient_upstream_cdna",
            "insufficient_downstream_cdna",
            "insufficient_cdna_on_both_sides",
        }
        for row in excluded_rows
    )

    print()
    print("Run complete")
    print(f"  APAlyzer rows read: {len(apa_results)}")
    print(f"  Significant UP PAS rows: {len(qualifying_rows)}")
    print(
        "  Unique significant UP genes: "
        f"{qualifying_rows['_gene_symbol_clean'].nunique()}"
    )
    print(f"  Intronic dIPAs retained: {intronic_count}")
    print(f"  Exonic dIPAs retained after exon removal: {exonic_count}")
    if legacy_strand_mismatch_count:
        print(
            "  WARNING: Legacy PASid strand disagreed with the current "
            f"canonical transcript for {legacy_strand_mismatch_count} gene(s)."
        )
    if legacy_strand_unavailable_count:
        print(
            "  Legacy PASid strand was unavailable for "
            f"{legacy_strand_unavailable_count} gene(s)."
        )
    print(f"  Genes excluded for short cDNA sides: {short_cdna_exclusions}")
    print(f"  Final genes used in every profile: {len(gene_models)}")
    print(
        "  Treatments plotted: "
        + ", ".join(treatment_condition_order)
    )
    print(f"  Threads used: {thread_count}")
    print(f"  Results written to: {output_directory}")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
