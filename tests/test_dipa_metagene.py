import json
import os
import subprocess
import sys
from pathlib import Path

import gffutils
import matplotlib.image as mpimg
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
            'chr1\ttest\tgene\t2601\t2990\t.\t+\t.\t'
            'gene_id "G_SHARED_1"; gene_name "SharedGene";'
        ),
        (
            'chr1\ttest\ttranscript\t2601\t2990\t.\t+\t.\t'
            'gene_id "G_SHARED_1"; transcript_id "TX_SHARED_1"; '
            'gene_name "SharedGene"; tag "Ensembl_canonical";'
        ),
        (
            'chr1\ttest\texon\t2601\t2700\t.\t+\t.\t'
            'gene_id "G_SHARED_1"; transcript_id "TX_SHARED_1"; '
            'gene_name "SharedGene"; exon_id "E_SHARED_1A";'
        ),
        (
            'chr1\ttest\texon\t2751\t2850\t.\t+\t.\t'
            'gene_id "G_SHARED_1"; transcript_id "TX_SHARED_1"; '
            'gene_name "SharedGene"; exon_id "E_SHARED_1B";'
        ),
        (
            'chr1\ttest\texon\t2891\t2990\t.\t+\t.\t'
            'gene_id "G_SHARED_1"; transcript_id "TX_SHARED_1"; '
            'gene_name "SharedGene"; exon_id "E_SHARED_1C";'
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
            'chr1\ttest\tgene\t3601\t3990\t.\t+\t.\t'
            'gene_id "G_SHARED_2"; gene_name "SharedGene";'
        ),
        (
            'chr1\ttest\ttranscript\t3601\t3990\t.\t+\t.\t'
            'gene_id "G_SHARED_2"; transcript_id "TX_SHARED_2"; '
            'gene_name "SharedGene"; tag "Ensembl_canonical";'
        ),
        (
            'chr1\ttest\texon\t3601\t3700\t.\t+\t.\t'
            'gene_id "G_SHARED_2"; transcript_id "TX_SHARED_2"; '
            'gene_name "SharedGene"; exon_id "E_SHARED_2A";'
        ),
        (
            'chr1\ttest\texon\t3751\t3850\t.\t+\t.\t'
            'gene_id "G_SHARED_2"; transcript_id "TX_SHARED_2"; '
            'gene_name "SharedGene"; exon_id "E_SHARED_2B";'
        ),
        (
            'chr1\ttest\texon\t3891\t3990\t.\t+\t.\t'
            'gene_id "G_SHARED_2"; transcript_id "TX_SHARED_2"; '
            'gene_name "SharedGene"; exon_id "E_SHARED_2C";'
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


def write_infinite_bigwig(bigwig_path, chromosome_length=5000):
    bigwig_file = pyBigWig.open(str(bigwig_path), "w")
    bigwig_file.addHeader([("chr1", chromosome_length)])
    bigwig_file.addEntries(
        ["chr1"],
        [0],
        ends=[chromosome_length],
        values=[float("inf")],
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
            "gene_symbol": "SharedGene",
            "PASid": "chr1:+:shared_coordinate_1",
            "RED": 3.5,
            "pvalue": 0.001,
            "p_adj": 0.01,
            "APAreg": "UP",
            "APAreg_original": "UP",
            "significance_for_plot": 0.01,
            "chr": "chr1",
            "start": 2869,
            "end": 2870,
            "browser_start_1based": 2870,
            "browser_end_1based": 2870,
            "browser_region": "chr1:2870-2870",
            "pas_coordinate_match": "TRUE",
        },
        {
            "gene_symbol": "SharedGene",
            "PASid": "chr1:+:shared_coordinate_2",
            "RED": 3.0,
            "pvalue": 0.001,
            "p_adj": 0.01,
            "APAreg": "UP",
            "APAreg_original": "UP",
            "significance_for_plot": 0.01,
            "chr": "chr1",
            "start": 3869,
            "end": 3870,
            "browser_start_1based": 3870,
            "browser_end_1based": 3870,
            "browser_region": "chr1:3870-3870",
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


def make_basic_paired_sample_rows(tmp_path):
    write_constant_bigwig(tmp_path / "basic_control.bw", 1)
    write_constant_bigwig(tmp_path / "basic_treatment.bw", 2)
    return [
        {
            "sample_id": "control_1",
            "condition": "control",
            "replicate": "1",
            "role": "control",
            "pair_id": "1",
            "bigwig": "basic_control.bw",
        },
        {
            "sample_id": "drugA_1",
            "condition": "drugA",
            "replicate": "1",
            "role": "treatment",
            "pair_id": "1",
            "bigwig": "basic_treatment.bw",
        },
    ]


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


def test_gtf_database_rejects_a_different_source_gtf(tmp_path):
    gtf_path = tmp_path / "annotation.gtf"
    database_path = tmp_path / "annotation.db"
    write_synthetic_gtf(gtf_path)

    database = dipa_metagene.create_or_open_gtf_database(
        gtf_path,
        database_path,
        False,
    )
    database.conn.close()

    metadata_path = tmp_path / "annotation.db.metadata.json"
    assert metadata_path.is_file()

    gtf_path.write_text(
        gtf_path.read_text(encoding="utf-8") + "# changed annotation\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="different GTF"):
        dipa_metagene.create_or_open_gtf_database(
            gtf_path,
            database_path,
            False,
        )


def test_gtf_database_rejects_missing_source_metadata(tmp_path):
    gtf_path = tmp_path / "annotation.gtf"
    database_path = tmp_path / "annotation.db"
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

    with pytest.raises(ValueError, match="no source metadata"):
        dipa_metagene.create_or_open_gtf_database(
            gtf_path,
            database_path,
            False,
        )


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
        "_bigwig_mode": "combined",
        "_combined_bigwig_path": str(bigwig_path),
    }

    _, profiles = dipa_metagene.extract_sample_profiles(
        sample_record,
        [minus_model],
        2,
    )

    assert profiles.shape == (1, 4)
    assert profiles[0, 0] > profiles[0, 1]
    assert profiles[0, 2] > profiles[0, 3]


def test_strand_specific_bigwigs_follow_annotation_strand(tmp_path):
    _, database = open_synthetic_gtf_database(tmp_path)

    plus_row = make_apa_series("GenePlus", "chr1:+:old", 350)
    plus_model, reason, _ = dipa_metagene.prepare_cdna_model(
        plus_row,
        database["TX_PLUS"],
        database,
        2,
    )
    assert reason is None

    minus_row = make_apa_series("GeneMinus", "chr1:-:old", 2350)
    minus_model, reason, _ = dipa_metagene.prepare_cdna_model(
        minus_row,
        database["TX_MINUS"],
        database,
        2,
    )
    assert reason is None

    plus_bigwig_path = tmp_path / "plus.bw"
    minus_bigwig_path = tmp_path / "minus.bw"
    write_constant_bigwig(plus_bigwig_path, 2)
    write_constant_bigwig(minus_bigwig_path, 7)

    sample_record = {
        "sample_id": "stranded_sample",
        "_bigwig_mode": "strand_specific",
        "_plus_bigwig_path": str(plus_bigwig_path),
        "_minus_bigwig_path": str(minus_bigwig_path),
    }

    _, profiles = dipa_metagene.extract_sample_profiles(
        sample_record,
        [plus_model, minus_model],
        2,
    )

    assert np.allclose(profiles[0], 2)
    assert np.allclose(profiles[1], 7)


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
        "_bigwig_mode": "combined",
        "_combined_bigwig_path": str(bigwig_path),
    }

    with pytest.raises(ValueError, match="Negative BigWig value"):
        dipa_metagene.extract_sample_profiles(
            sample_record,
            [plus_model],
            2,
        )


def test_infinite_bigwig_signal_is_rejected(tmp_path):
    _, database = open_synthetic_gtf_database(tmp_path)
    plus_row = make_apa_series("GenePlus", "chr1:+:old", 350)
    plus_model, reason, _ = dipa_metagene.prepare_cdna_model(
        plus_row,
        database["TX_PLUS"],
        database,
        2,
    )
    assert reason is None

    bigwig_path = tmp_path / "infinite.bw"
    write_infinite_bigwig(bigwig_path)
    sample_record = {
        "sample_id": "infinite_sample",
        "_bigwig_mode": "combined",
        "_combined_bigwig_path": str(bigwig_path),
    }

    with pytest.raises(ValueError, match="Infinite BigWig value"):
        dipa_metagene.extract_sample_profiles(
            sample_record,
            [plus_model],
            2,
        )


def test_blank_and_duplicate_replicates_are_rejected(tmp_path):
    blank_rows = make_basic_paired_sample_rows(tmp_path)
    blank_rows[1]["replicate"] = ""

    blank_run, _ = run_script(
        tmp_path,
        blank_rows,
        "blank_replicate_output",
    )
    assert blank_run.returncode == 1
    assert "blank replicate" in blank_run.stderr

    write_constant_bigwig(tmp_path / "duplicate_treatment.bw", 3)
    duplicate_rows = make_basic_paired_sample_rows(tmp_path)
    duplicate_rows.append(
        {
            "sample_id": "drugA_duplicate",
            "condition": "drugA",
            "replicate": "1",
            "role": "treatment",
            "pair_id": "1",
            "bigwig": "duplicate_treatment.bw",
        }
    )

    duplicate_run, _ = run_script(
        tmp_path,
        duplicate_rows,
        "duplicate_replicate_output",
    )
    assert duplicate_run.returncode == 1
    assert "Replicate labels must be unique" in duplicate_run.stderr


def test_invalid_bigwig_column_combinations_are_rejected(tmp_path):
    plus_only_rows = make_basic_paired_sample_rows(tmp_path)
    write_constant_bigwig(tmp_path / "treatment_plus.bw", 2)
    plus_only_rows[1]["bigwig"] = ""
    plus_only_rows[1]["plus_bigwig"] = "treatment_plus.bw"

    plus_only_run, _ = run_script(
        tmp_path,
        plus_only_rows,
        "plus_only_output",
    )
    assert plus_only_run.returncode == 1
    assert (
        "either bigwig or both plus_bigwig and minus_bigwig"
        in plus_only_run.stderr
    )

    combined_and_stranded_rows = make_basic_paired_sample_rows(tmp_path)
    write_constant_bigwig(tmp_path / "treatment_minus.bw", 2)
    combined_and_stranded_rows[1]["plus_bigwig"] = "treatment_plus.bw"
    combined_and_stranded_rows[1]["minus_bigwig"] = "treatment_minus.bw"

    combined_and_stranded_run, _ = run_script(
        tmp_path,
        combined_and_stranded_rows,
        "combined_and_stranded_output",
    )
    assert combined_and_stranded_run.returncode == 1
    assert (
        "both a combined BigWig and a strand-specific BigWig"
        in combined_and_stranded_run.stderr
    )


@pytest.mark.parametrize("invalid_pseudocount", ["nan", "inf"])
def test_nonfinite_pseudocount_is_rejected(
    tmp_path,
    invalid_pseudocount,
):
    sample_rows = make_basic_paired_sample_rows(tmp_path)
    completed, _ = run_script(
        tmp_path,
        sample_rows,
        f"invalid_pseudocount_{invalid_pseudocount}",
        extra_arguments=["--pseudocount", invalid_pseudocount],
    )

    assert completed.returncode == 1
    assert "pseudocount must be a finite number" in completed.stderr


@pytest.mark.parametrize(
    ("column_name", "invalid_value", "expected_message"),
    [
        ("RED", "inf", "contains an infinite or missing numeric value"),
        ("p_adj", "-0.1", "must contain values between 0 and 1"),
        ("p_adj", "1.1", "must contain values between 0 and 1"),
        (
            "p_adj",
            "inf",
            "contains an infinite numeric value",
        ),
        ("p_adj", "not_a_number", "contains a nonnumeric value"),
    ],
)
def test_invalid_apalyzer_numeric_values_are_rejected(
    tmp_path,
    column_name,
    invalid_value,
    expected_message,
):
    apa_path = tmp_path / "apa.tsv"
    write_apa_table(apa_path)
    apa_table = pd.read_csv(apa_path, sep="\t", dtype=str)
    apa_table.loc[0, column_name] = invalid_value
    apa_table.to_csv(apa_path, sep="\t", index=False)

    sample_rows = make_basic_paired_sample_rows(tmp_path)
    completed, _ = run_script(
        tmp_path,
        sample_rows,
        f"invalid_{column_name}",
    )

    assert completed.returncode == 1
    assert expected_message in completed.stderr


def test_missing_adjusted_pvalue_is_excluded(tmp_path):
    apa_path = tmp_path / "apa.tsv"
    write_apa_table(apa_path)
    apa_table = pd.read_csv(apa_path, sep="\t", dtype=str)
    apa_table.loc[0, "p_adj"] = "NA"
    apa_table.to_csv(apa_path, sep="\t", index=False)

    sample_rows = make_basic_paired_sample_rows(tmp_path)
    completed, output_path = run_script(
        tmp_path,
        sample_rows,
        "missing_adjusted_pvalue_output",
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout

    excluded = pd.read_csv(
        output_path / "excluded_genes.tsv",
        sep="\t",
        keep_default_na=False,
    )
    missing_rows = excluded[
        excluded["exclusion_reason"] == "missing_adjusted_pvalue"
    ]
    assert len(missing_rows) == 1
    assert missing_rows.iloc[0]["p_adj"] == "NA"

    selected = pd.read_csv(
        output_path / "selected_up_genes.tsv",
        sep="\t",
    )
    selected_plus = selected[
        selected["gene_symbol"] == "GenePlus"
    ].iloc[0]
    assert selected_plus["RED"] == pytest.approx(2.0)


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
        "synthetic.gtf.gffutils.db.metadata.json",
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
        "SharedGene",
    }
    assert len(selected) == 5
    shared_gene_ids = set(
        selected.loc[
            selected["gene_symbol"] == "SharedGene",
            "gene_id",
        ]
    )
    assert shared_gene_ids == {"G_SHARED_1", "G_SHARED_2"}
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

    per_gene = pd.read_csv(output_path / "per_gene_log2fc.tsv", sep="\t")
    replicate_profiles = pd.read_csv(
        output_path / "replicate_profiles.tsv",
        sep="\t",
    )
    assert set(per_gene["replicate"].astype(str)) == {"1", "2"}
    assert set(
        per_gene.loc[
            per_gene["gene_symbol"] == "SharedGene",
            "gene_id",
        ]
    ) == {"G_SHARED_1", "G_SHARED_2"}
    assert set(replicate_profiles["replicate"].astype(str)) == {"1", "2"}

    run_parameters = json.loads(
        (output_path / "run_parameters.json").read_text(encoding="utf-8")
    )
    assert run_parameters["parameters"]["paired_controls"] is True
    assert run_parameters["parameters"]["threads_used"] == 2
    assert run_parameters["filtering_counts"]["final_genes"] == 5
    assert (
        run_parameters["filtering_counts"]["unique_significant_up_genes"]
        == 5
    )
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


def test_full_run_accepts_mixed_bigwig_input_modes(tmp_path):
    write_constant_bigwig(tmp_path / "mixed_control.bw", 1)
    write_constant_bigwig(tmp_path / "mixed_treatment_plus.bw", 3)
    write_constant_bigwig(tmp_path / "mixed_treatment_minus.bw", 3)

    sample_rows = [
        {
            "sample_id": "control_1",
            "condition": "control",
            "replicate": "1",
            "role": "control",
            "pair_id": "1",
            "bigwig": "mixed_control.bw",
            "plus_bigwig": "",
            "minus_bigwig": "",
        },
        {
            "sample_id": "drugA_1",
            "condition": "drugA",
            "replicate": "1",
            "role": "treatment",
            "pair_id": "1",
            "bigwig": "",
            "plus_bigwig": "mixed_treatment_plus.bw",
            "minus_bigwig": "mixed_treatment_minus.bw",
        },
    ]

    completed, output_path = run_script(
        tmp_path,
        sample_rows,
        "mixed_bigwig_output",
    )

    assert completed.returncode == 0, completed.stderr + completed.stdout

    summary = pd.read_csv(
        output_path / "treatment_summary.tsv",
        sep="\t",
    )
    assert np.allclose(summary["mean_log2fc"], 1.0)

    run_parameters = json.loads(
        (output_path / "run_parameters.json").read_text(encoding="utf-8")
    )
    assert run_parameters["parameters"]["combined_bigwig_samples"] == 1
    assert (
        run_parameters["parameters"]["strand_specific_bigwig_samples"]
        == 1
    )


def test_plot_footer_fits_with_multiple_treatments(tmp_path):
    write_constant_bigwig(tmp_path / "footer_control.bw", 1)

    sample_rows = [
        {
            "sample_id": "control_1",
            "condition": "control",
            "replicate": "1",
            "role": "control",
            "pair_id": "1",
            "bigwig": "footer_control.bw",
        }
    ]

    treatment_conditions = [
        "AS_PHA_4hr",
        "AS_PHA_8hr",
        "AS_PHA_24hr",
        "AS_PHA_48hr",
    ]
    for condition in treatment_conditions:
        bigwig_name = f"{condition}.bw"
        write_constant_bigwig(tmp_path / bigwig_name, 3)
        sample_rows.append(
            {
                "sample_id": f"{condition}_1",
                "condition": condition,
                "replicate": "1",
                "role": "treatment",
                "pair_id": "1",
                "bigwig": bigwig_name,
            }
        )

    completed, output_path = run_script(
        tmp_path,
        sample_rows,
        "footer_output",
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout

    plot_image = mpimg.imread(output_path / "combined_metagene.png")
    bottom_region = plot_image[
        int(plot_image.shape[0] * 0.93) :,
        :,
        :3,
    ]
    text_pixels = np.any(bottom_region < 0.85, axis=2)
    horizontal_positions = np.where(text_pixels)[1]

    assert len(horizontal_positions) > 0
    assert horizontal_positions.min() > 5
    assert horizontal_positions.max() < plot_image.shape[1] - 6
