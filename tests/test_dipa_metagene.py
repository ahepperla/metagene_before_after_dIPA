import json
import os
import subprocess
import sys
from pathlib import Path

import gffutils
import numpy as np
import pandas as pd
import pyBigWig
import pytest


PROJECT_DIRECTORY = Path(__file__).resolve().parents[1]
SCRIPT_PATH = PROJECT_DIRECTORY / "dipa_metagene.py"
sys.path.insert(0, str(PROJECT_DIRECTORY))

import dipa_metagene


def write_synthetic_gtf(gtf_path):
    """Write small plus- and minus-strand canonical transcripts."""

    gtf_lines = [
        (
            'chr1\ttest\tgene\t1\t520\t.\t+\t.\t'
            'gene_id "G_PLUS"; gene_name "GenePlus";'
        ),
        (
            'chr1\ttest\ttranscript\t1\t520\t.\t+\t.\t'
            'gene_id "G_PLUS"; transcript_id "TX_PLUS"; '
            'gene_name "GenePlus"; tag "Ensembl_canonical";'
        ),
        (
            'chr1\ttest\texon\t1\t120\t.\t+\t.\t'
            'gene_id "G_PLUS"; transcript_id "TX_PLUS"; '
            'gene_name "GenePlus"; exon_id "E_PLUS_1";'
        ),
        (
            'chr1\ttest\texon\t201\t320\t.\t+\t.\t'
            'gene_id "G_PLUS"; transcript_id "TX_PLUS"; '
            'gene_name "GenePlus"; exon_id "E_PLUS_2";'
        ),
        (
            'chr1\ttest\texon\t401\t520\t.\t+\t.\t'
            'gene_id "G_PLUS"; transcript_id "TX_PLUS"; '
            'gene_name "GenePlus"; exon_id "E_PLUS_3";'
        ),
        (
            'chr1\ttest\tgene\t1001\t1520\t.\t+\t.\t'
            'gene_id "G_EXONIC"; gene_name "GeneExonic";'
        ),
        (
            'chr1\ttest\ttranscript\t1001\t1520\t.\t+\t.\t'
            'gene_id "G_EXONIC"; transcript_id "TX_EXONIC"; '
            'gene_name "GeneExonic"; tag "Ensembl_canonical";'
        ),
        (
            'chr1\ttest\texon\t1001\t1120\t.\t+\t.\t'
            'gene_id "G_EXONIC"; transcript_id "TX_EXONIC"; '
            'gene_name "GeneExonic"; exon_id "E_EXONIC_1";'
        ),
        (
            'chr1\ttest\texon\t1201\t1320\t.\t+\t.\t'
            'gene_id "G_EXONIC"; transcript_id "TX_EXONIC"; '
            'gene_name "GeneExonic"; exon_id "E_EXONIC_2";'
        ),
        (
            'chr1\ttest\texon\t1401\t1520\t.\t+\t.\t'
            'gene_id "G_EXONIC"; transcript_id "TX_EXONIC"; '
            'gene_name "GeneExonic"; exon_id "E_EXONIC_3";'
        ),
        (
            'chr1\ttest\tgene\t2001\t2520\t.\t-\t.\t'
            'gene_id "G_MINUS"; gene_name "GeneMinus";'
        ),
        (
            'chr1\ttest\ttranscript\t2001\t2520\t.\t-\t.\t'
            'gene_id "G_MINUS"; transcript_id "TX_MINUS"; '
            'gene_name "GeneMinus"; tag "Ensembl_canonical";'
        ),
        (
            'chr1\ttest\texon\t2001\t2120\t.\t-\t.\t'
            'gene_id "G_MINUS"; transcript_id "TX_MINUS"; '
            'gene_name "GeneMinus"; exon_id "E_MINUS_1";'
        ),
        (
            'chr1\ttest\texon\t2201\t2320\t.\t-\t.\t'
            'gene_id "G_MINUS"; transcript_id "TX_MINUS"; '
            'gene_name "GeneMinus"; exon_id "E_MINUS_2";'
        ),
        (
            'chr1\ttest\texon\t2401\t2520\t.\t-\t.\t'
            'gene_id "G_MINUS"; transcript_id "TX_MINUS"; '
            'gene_name "GeneMinus"; exon_id "E_MINUS_3";'
        ),
        (
            'chr1\ttest\tgene\t3001\t3500\t.\t+\t.\t'
            'gene_id "G_SHORT"; gene_name "GeneShort";'
        ),
        (
            'chr1\ttest\ttranscript\t3001\t3500\t.\t+\t.\t'
            'gene_id "G_SHORT"; transcript_id "TX_SHORT"; '
            'gene_name "GeneShort"; tag "Ensembl_canonical";'
        ),
        (
            'chr1\ttest\texon\t3001\t3099\t.\t+\t.\t'
            'gene_id "G_SHORT"; transcript_id "TX_SHORT"; '
            'gene_name "GeneShort"; exon_id "E_SHORT_1";'
        ),
        (
            'chr1\ttest\texon\t3201\t3300\t.\t+\t.\t'
            'gene_id "G_SHORT"; transcript_id "TX_SHORT"; '
            'gene_name "GeneShort"; exon_id "E_SHORT_2";'
        ),
        (
            'chr1\ttest\texon\t3401\t3500\t.\t+\t.\t'
            'gene_id "G_SHORT"; transcript_id "TX_SHORT"; '
            'gene_name "GeneShort"; exon_id "E_SHORT_3";'
        ),
        (
            'chr1\ttest\tgene\t4001\t4520\t.\t+\t.\t'
            'gene_id "G_NO_CANONICAL"; gene_name "GeneNoCanonical";'
        ),
        (
            'chr1\ttest\ttranscript\t4001\t4520\t.\t+\t.\t'
            'gene_id "G_NO_CANONICAL"; transcript_id "TX_NO_CANONICAL"; '
            'gene_name "GeneNoCanonical";'
        ),
        (
            'chr1\ttest\texon\t4001\t4120\t.\t+\t.\t'
            'gene_id "G_NO_CANONICAL"; transcript_id "TX_NO_CANONICAL"; '
            'gene_name "GeneNoCanonical"; exon_id "E_NO_CANONICAL_1";'
        ),
        (
            'chr1\ttest\texon\t4401\t4520\t.\t+\t.\t'
            'gene_id "G_NO_CANONICAL"; transcript_id "TX_NO_CANONICAL"; '
            'gene_name "GeneNoCanonical"; exon_id "E_NO_CANONICAL_2";'
        ),
    ]
    gtf_path.write_text("\n".join(gtf_lines) + "\n", encoding="utf-8")


def open_synthetic_gtf_database(tmp_path):
    gtf_path = tmp_path / "synthetic.gtf"
    database_path = tmp_path / "synthetic.db"
    write_synthetic_gtf(gtf_path)
    gffutils.create_db(
        str(gtf_path),
        dbfn=str(database_path),
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
        disable_infer_genes=True,
        disable_infer_transcripts=True,
    )
    return gtf_path, gffutils.FeatureDB(str(database_path), keep_order=True)


def make_apa_series(gene_symbol, pas_id, dipa_position, red=5.0):
    """Create an APA row after numeric and coordinate validation."""

    return pd.Series(
        {
            "gene_symbol": gene_symbol,
            "PASid": pas_id,
            "RED": str(red),
            "p_adj": "0.01",
            "APAreg": "UP",
            "APAreg_original": "UP",
            "chr": "chr1",
            "start": str(dipa_position - 1),
            "end": str(dipa_position),
            "browser_start_1based": str(dipa_position),
            "pas_coordinate_match": "TRUE",
            "_RED_numeric": red,
            "_p_adj_numeric": 0.01,
            "_start_numeric": dipa_position - 1,
            "_browser_start_numeric": dipa_position,
            "_legacy_pas_strand": dipa_metagene.parse_pas_strand(pas_id),
        }
    )


def write_constant_bigwig(bigwig_path, value, chromosome_length=5000):
    bigwig_file = pyBigWig.open(str(bigwig_path), "w")
    bigwig_file.addHeader([("chr1", chromosome_length)])
    bigwig_file.addEntries(
        ["chr1"],
        [0],
        ends=[chromosome_length],
        values=[float(value)],
    )
    bigwig_file.close()


def write_position_bigwig(bigwig_path, chromosome_length=5000):
    """Write one value per base so strand reversal is directly testable."""

    starts = list(range(chromosome_length))
    ends = [start + 1 for start in starts]
    values = [float(start + 1) for start in starts]
    bigwig_file = pyBigWig.open(str(bigwig_path), "w")
    bigwig_file.addHeader([("chr1", chromosome_length)])
    bigwig_file.addEntries(
        ["chr1"] * chromosome_length,
        starts,
        ends=ends,
        values=values,
    )
    bigwig_file.close()


def write_negative_bigwig(bigwig_path, chromosome_length=5000):
    bigwig_file = pyBigWig.open(str(bigwig_path), "w")
    bigwig_file.addHeader([("chr1", chromosome_length)])
    bigwig_file.addEntries(
        ["chr1"],
        [0],
        ends=[chromosome_length],
        values=[-1.0],
    )
    bigwig_file.close()


def write_apa_table(apa_path):
    rows = [
        {
            "gene_symbol": "GenePlus",
            "PASid": "chr1:+:old_coordinate_1",
            "RED": 5.0,
            "pvalue": 0.001,
            "p_adj": 0.01,
            "APAreg": "UP",
            "APAreg_original": "UP",
            "significance_for_plot": 0.01,
            "chr": "chr1",
            "start": 349,
            "end": 350,
            "browser_start_1based": 350,
            "browser_end_1based": 350,
            "browser_region": "chr1:350-350",
            "pas_coordinate_match": "TRUE",
        },
        {
            "gene_symbol": "GenePlus",
            "PASid": "chr1:+:old_coordinate_2",
            "RED": 2.0,
            "pvalue": 0.001,
            "p_adj": 0.01,
            "APAreg": "UP",
            "APAreg_original": "UP",
            "significance_for_plot": 0.01,
            "chr": "chr1",
            "start": 359,
            "end": 360,
            "browser_start_1based": 360,
            "browser_end_1based": 360,
            "browser_region": "chr1:360-360",
            "pas_coordinate_match": "TRUE",
        },
        {
            "gene_symbol": "GeneExonic",
            "PASid": "legacy_exonic_label_without_strand",
            "RED": 4.0,
            "pvalue": 0.001,
            "p_adj": 0.01,
            "APAreg": "UP",
            "APAreg_original": "UP",
            "significance_for_plot": 0.01,
            "chr": "chr1",
            "start": 1249,
            "end": 1250,
            "browser_start_1based": 1250,
            "browser_end_1based": 1250,
            "browser_region": "chr1:1250-1250",
            "pas_coordinate_match": "TRUE",
        },
        {
            "gene_symbol": "GeneMinus",
            # This deliberately disagrees with the current minus-strand GTF.
            # The legacy PASid must remain QC-only after coordinate liftOver.
            "PASid": "chr1:+:old_coordinate_4",
            "RED": 3.0,
            "pvalue": 0.001,
            "p_adj": 0.01,
            "APAreg": "UP",
            "APAreg_original": "UP",
            "significance_for_plot": 0.01,
            "chr": "chr1",
            "start": 2349,
            "end": 2350,
            "browser_start_1based": 2350,
            "browser_end_1based": 2350,
            "browser_region": "chr1:2350-2350",
            "pas_coordinate_match": "TRUE",
        },
        {
            "gene_symbol": "GeneNoCanonical",
            "PASid": "chr1:+:old_coordinate_5",
            "RED": 3.0,
            "pvalue": 0.001,
            "p_adj": 0.01,
            "APAreg": "UP",
            "APAreg_original": "UP",
            "significance_for_plot": 0.01,
            "chr": "chr1",
            "start": 4299,
            "end": 4300,
            "browser_start_1based": 4300,
            "browser_end_1based": 4300,
            "browser_region": "chr1:4300-4300",
            "pas_coordinate_match": "TRUE",
        },
        {
            "gene_symbol": "DownGene",
            "PASid": "chr1:+:old_coordinate_6",
            "RED": -2.0,
            "pvalue": 0.001,
            "p_adj": 0.01,
            "APAreg": "DN",
            "APAreg_original": "DN",
            "significance_for_plot": 0.01,
            "chr": "chr1",
            "start": 599,
            "end": 600,
            "browser_start_1based": 600,
            "browser_end_1based": 600,
            "browser_region": "chr1:600-600",
            "pas_coordinate_match": "TRUE",
        },
        {
            "gene_symbol": "NotSignificant",
            "PASid": "chr1:+:old_coordinate_7",
            "RED": 2.0,
            "pvalue": 0.5,
            "p_adj": 0.5,
            "APAreg": "UP",
            "APAreg_original": "UP",
            "significance_for_plot": 0.5,
            "chr": "chr1",
            "start": 699,
            "end": 700,
            "browser_start_1based": 700,
            "browser_end_1based": 700,
            "browser_region": "chr1:700-700",
            "pas_coordinate_match": "TRUE",
        },
    ]
    pd.DataFrame(rows).to_csv(apa_path, sep="\t", index=False)


def run_script(
    tmp_path,
    sample_rows,
    output_name,
    extra_arguments=None,
):
    gtf_path = tmp_path / "synthetic.gtf"
    apa_path = tmp_path / "apa.tsv"
    sample_path = tmp_path / f"{output_name}_samples.tsv"
    output_path = tmp_path / output_name

    if not gtf_path.exists():
        write_synthetic_gtf(gtf_path)
    if not apa_path.exists():
        write_apa_table(apa_path)

    pd.DataFrame(sample_rows).to_csv(sample_path, sep="\t", index=False)

    command = [
        sys.executable,
        str(SCRIPT_PATH),
        "--apa-results",
        str(apa_path),
        "--gtf",
        str(gtf_path),
        "--samples",
        str(sample_path),
        "--output-dir",
        str(output_path),
        "--pseudocount",
        "1",
        "--adjusted-pvalue",
        "0.05",
        "--bins-per-side",
        "10",
        "--min-side-bp",
        "10",
        "--threads",
        "2",
    ]
    if extra_arguments:
        command.extend(extra_arguments)

    completed = subprocess.run(
        command,
        cwd=PROJECT_DIRECTORY,
        text=True,
        capture_output=True,
        check=False,
    )
    return completed, output_path


def test_parse_pas_strand_and_gene_id_normalization():
    assert dipa_metagene.parse_pas_strand("chr1:+:123") == "+"
    assert dipa_metagene.parse_pas_strand("chrX:-:999") == "-"
    assert dipa_metagene.parse_pas_strand("bad_value") is None
    assert dipa_metagene.normalize_gene_id("ENSMUSG000001.12") == (
        "ENSMUSG000001"
    )


def test_bin_cdna_signal_uses_real_nonempty_bases():
    values = np.arange(200, dtype=float)
    binned = dipa_metagene.bin_cdna_signal(values, 100)

    assert len(binned) == 100
    assert binned[0] == pytest.approx(0.5)
    assert binned[-1] == pytest.approx(198.5)

    with pytest.raises(ValueError, match="nonempty bins"):
        dipa_metagene.bin_cdna_signal(np.arange(99), 100)


def test_thread_count_respects_request_and_slurm(monkeypatch):
    assert dipa_metagene.choose_thread_count(8, 3) == 3

    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.delenv("PBS_NP", raising=False)
    monkeypatch.delenv("NSLOTS", raising=False)
    assert dipa_metagene.choose_thread_count(0, 5) == 2

    with pytest.raises(ValueError, match="zero or a positive"):
        dipa_metagene.choose_thread_count(-1, 5)


def test_intronic_exonic_minus_and_short_cdna_models(tmp_path):
    _, database = open_synthetic_gtf_database(tmp_path)

    plus_row = make_apa_series("GenePlus", "chr1:+:old", 350)
    plus_model, reason, _ = dipa_metagene.prepare_cdna_model(
        plus_row,
        database["TX_PLUS"],
        database,
        100,
    )
    assert reason is None
    assert plus_model["dipa_context"] == "intronic"
    assert [
        exon["exon_id"] for exon in plus_model["upstream_exons"]
    ] == ["E_PLUS_1", "E_PLUS_2"]
    assert [
        exon["exon_id"] for exon in plus_model["downstream_exons"]
    ] == ["E_PLUS_3"]

    exonic_row = make_apa_series("GeneExonic", "chr1:+:old", 1250)
    exonic_model, reason, _ = dipa_metagene.prepare_cdna_model(
        exonic_row,
        database["TX_EXONIC"],
        database,
        100,
    )
    assert reason is None
    assert exonic_model["dipa_context"] == "exonic_exon_removed"
    assert exonic_model["removed_exon_id"] == "E_EXONIC_2"
    assert [
        exon["exon_id"] for exon in exonic_model["upstream_exons"]
    ] == ["E_EXONIC_1"]
    assert [
        exon["exon_id"] for exon in exonic_model["downstream_exons"]
    ] == ["E_EXONIC_3"]

    minus_row = make_apa_series("GeneMinus", "chr1:-:old", 2350)
    minus_model, reason, _ = dipa_metagene.prepare_cdna_model(
        minus_row,
        database["TX_MINUS"],
        database,
        100,
    )
    assert reason is None
    assert [
        exon["exon_id"] for exon in minus_model["upstream_exons"]
    ] == ["E_MINUS_3"]
    assert [
        exon["exon_id"] for exon in minus_model["downstream_exons"]
    ] == ["E_MINUS_2", "E_MINUS_1"]

    short_row = make_apa_series("GeneShort", "chr1:+:old", 3150)
    short_model, reason, details = dipa_metagene.prepare_cdna_model(
        short_row,
        database["TX_SHORT"],
        database,
        100,
    )
    assert short_model is None
    assert reason == "insufficient_upstream_cdna"
    assert "99 bp" in details


def test_minus_strand_bigwig_signal_is_reversed(tmp_path):
    _, database = open_synthetic_gtf_database(tmp_path)
    minus_row = make_apa_series("GeneMinus", "chr1:-:old", 2350)
    minus_model, reason, _ = dipa_metagene.prepare_cdna_model(
        minus_row,
        database["TX_MINUS"],
        database,
        2,
    )
    assert reason is None

    bigwig_path = tmp_path / "positions.bw"
    write_position_bigwig(bigwig_path)
    sample_record = {
        "sample_id": "sample",
        "_bigwig_path": str(bigwig_path),
    }

    _, profiles = dipa_metagene.extract_sample_profiles(
        sample_record,
        [minus_model],
        2,
    )

    assert profiles.shape == (1, 4)
    assert profiles[0, 0] > profiles[0, 1]
    assert profiles[0, 2] > profiles[0, 3]


def test_negative_bigwig_signal_is_rejected(tmp_path):
    _, database = open_synthetic_gtf_database(tmp_path)
    plus_row = make_apa_series("GenePlus", "chr1:+:old", 350)
    plus_model, reason, _ = dipa_metagene.prepare_cdna_model(
        plus_row,
        database["TX_PLUS"],
        database,
        2,
    )
    assert reason is None

    bigwig_path = tmp_path / "negative.bw"
    write_negative_bigwig(bigwig_path)
    sample_record = {
        "sample_id": "negative_sample",
        "_bigwig_path": str(bigwig_path),
    }

    with pytest.raises(ValueError, match="Negative BigWig value"):
        dipa_metagene.extract_sample_profiles(
            sample_record,
            [plus_model],
            2,
        )


def test_full_paired_run_with_multiple_treatments(tmp_path):
    signal_values = {
        "control_1": 1,
        "control_2": 2,
        "drugA_1": 3,
        "drugA_2": 5,
        "drugB_1": 1,
        "drugB_2": 2,
    }
    sample_rows = []

    for sample_id, signal_value in signal_values.items():
        bigwig_name = f"{sample_id}.bw"
        write_constant_bigwig(tmp_path / bigwig_name, signal_value)

        if sample_id.startswith("control"):
            role = "control"
            condition = "control"
        else:
            role = "treatment"
            condition = sample_id.split("_")[0]

        replicate = sample_id.rsplit("_", 1)[1]
        sample_rows.append(
            {
                "sample_id": sample_id,
                "condition": condition,
                "replicate": replicate,
                "role": role,
                "pair_id": replicate,
                "bigwig": bigwig_name,
            }
        )

    completed, output_path = run_script(
        tmp_path,
        sample_rows,
        "paired_output",
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout

    expected_files = [
        "combined_metagene.png",
        "combined_metagene.pdf",
        "selected_up_genes.tsv",
        "excluded_genes.tsv",
        "per_gene_log2fc.tsv",
        "replicate_profiles.tsv",
        "treatment_summary.tsv",
        "run_parameters.json",
    ]
    for filename in expected_files:
        assert (output_path / filename).is_file()
        assert (output_path / filename).stat().st_size > 0

    selected = pd.read_csv(
        output_path / "selected_up_genes.tsv",
        sep="\t",
    )
    assert set(selected["gene_symbol"]) == {
        "GenePlus",
        "GeneExonic",
        "GeneMinus",
    }
    selected_plus = selected[selected["gene_symbol"] == "GenePlus"].iloc[0]
    assert selected_plus["RED"] == pytest.approx(5.0)
    selected_exonic = selected[
        selected["gene_symbol"] == "GeneExonic"
    ].iloc[0]
    assert selected_exonic["dipa_context"] == "exonic_exon_removed"
    assert selected_exonic["removed_exon_id"] == "E_EXONIC_2"
    assert pd.isna(selected_exonic["legacy_pas_strand"])
    assert (
        selected_exonic["legacy_strand_matches_annotation"]
        == "unavailable"
    )

    selected_minus = selected[
        selected["gene_symbol"] == "GeneMinus"
    ].iloc[0]
    assert selected_minus["strand"] == "-"
    assert selected_minus["strand_source"] == (
        "Ensembl_canonical_transcript"
    )
    assert selected_minus["legacy_pas_strand"] == "+"
    assert selected_minus["legacy_strand_matches_annotation"] == "no"

    excluded = pd.read_csv(output_path / "excluded_genes.tsv", sep="\t")
    exclusion_reasons = set(excluded["exclusion_reason"])
    assert "additional_significant_PAS_for_gene" in exclusion_reasons
    assert "not_UP" in exclusion_reasons
    assert "not_significant" in exclusion_reasons
    assert "no_matching_Ensembl_canonical_transcript" in exclusion_reasons

    summary = pd.read_csv(output_path / "treatment_summary.tsv", sep="\t")
    drug_a = summary[summary["condition"] == "drugA"]
    drug_b = summary[summary["condition"] == "drugB"]
    assert np.allclose(drug_a["mean_log2fc"], 1.0)
    assert np.allclose(drug_a["sd_log2fc"], 0.0)
    assert np.allclose(drug_b["mean_log2fc"], 0.0)
    assert np.allclose(drug_b["sd_log2fc"], 0.0)

    run_parameters = json.loads(
        (output_path / "run_parameters.json").read_text(encoding="utf-8")
    )
    assert run_parameters["parameters"]["paired_controls"] is True
    assert run_parameters["parameters"]["threads_used"] == 2
    assert run_parameters["filtering_counts"]["final_genes"] == 3
    assert (
        run_parameters["filtering_counts"][
            "legacy_pas_strand_mismatches"
        ]
        == 1
    )
    assert (
        run_parameters["filtering_counts"][
            "legacy_pas_strand_unavailable"
        ]
        == 1
    )
    assert "Legacy PASid strand disagreed" in completed.stdout


def test_full_unpaired_run_uses_mean_control(tmp_path):
    signal_values = {
        "control_1": 1,
        "control_2": 3,
        "drugA_1": 5,
    }
    sample_rows = []

    for sample_id, signal_value in signal_values.items():
        bigwig_name = f"unpaired_{sample_id}.bw"
        write_constant_bigwig(tmp_path / bigwig_name, signal_value)

        role = "control" if sample_id.startswith("control") else "treatment"
        condition = "control" if role == "control" else "drugA"
        sample_rows.append(
            {
                "sample_id": sample_id,
                "condition": condition,
                "replicate": sample_id.rsplit("_", 1)[1],
                "role": role,
                "pair_id": "",
                "bigwig": bigwig_name,
            }
        )

    completed, output_path = run_script(
        tmp_path,
        sample_rows,
        "unpaired_output",
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert "Controls are unpaired" in completed.stdout

    summary = pd.read_csv(output_path / "treatment_summary.tsv", sep="\t")
    assert np.allclose(summary["mean_log2fc"], 1.0)
    assert summary["sd_log2fc"].isna().all()

    run_parameters = json.loads(
        (output_path / "run_parameters.json").read_text(encoding="utf-8")
    )
    assert run_parameters["parameters"]["paired_controls"] is False
