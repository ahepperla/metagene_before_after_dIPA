# dIPA Metagene Plot

`dipa_metagene.py` creates an exon-only metagene plot around differential
polyadenylation sites reported by APAlyzer.

The script:

1. Keeps significant APAlyzer events in the `UP` direction.
2. Chooses one dIPA per gene.
3. Uses the transcript tagged `Ensembl_canonical` in a GENCODE GTF.
4. Builds cDNA profiles from exons rather than the full genomic gene body.
5. Reads normalized signal from per-replicate BigWigs.
6. Calculates treatment-versus-control log2 fold-change.
7. Draws one line per treatment with replicate standard deviation shading.

The code is intentionally written as a readable script. It uses descriptive
variables, explicit validation, and only a small number of helpers for work
that repeats across genes or samples.

## Important Interpretation

### RED is used for event selection

APAlyzer's `RED` column is Relative Expression Difference. Positive values
indicate increased relative IPA usage in the treatment comparison.

This script uses RED to:

- confirm that an `UP` row has a positive direction
- choose the highest-RED dIPA when a gene has several significant UP sites

RED is not plotted. The y-axis is calculated independently from the BigWigs.

### The dIPA is a conceptual center

For an intronic dIPA, the script places exons before the containing intron on
the left and exons after the intron on the right. Intronic sequence is not
included.

For an exonic dIPA, the entire containing exon is removed. Exons before it are
placed on the left and exons after it are placed on the right.

The x-axis therefore represents:

```text
0% cDNA  --------  dIPA  --------  100% cDNA
          upstream        downstream
```

It does not display signal from the intronic PAS nucleotide or from an exon
that was removed because it contained the dIPA.

## Installation

Python 3.10 through 3.12 is recommended for broad HPC package compatibility.

### Conda or Mamba

This is usually the simplest HPC installation:

```bash
mamba env create -f environment.yml
conda activate dipa-metagene
```

Use `conda env create` if Mamba is unavailable.

### Python virtual environment

```bash
module load python/3.11
python -m venv dipa-metagene-env
source dipa-metagene-env/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The exact `module load` command depends on the HPC.

## Input Files

All three input files and all BigWigs must use the same genome assembly.
Chromosome names must also agree, such as `chr1` everywhere rather than `chr1`
in one file and `1` in another.

### APAlyzer table

The APAlyzer input is a tab-separated table. Required columns are:

```text
gene_symbol
PASid
RED
p_adj
APAreg
APAreg_original
chr
start
end
browser_start_1based
pas_coordinate_match
```

The script keeps rows satisfying:

```text
p_adj <= --adjusted-pvalue
APAreg == UP
APAreg_original == UP
RED > 0
pas_coordinate_match == TRUE
```

It uses the largest positive RED for genes with multiple qualifying sites.
Ties are resolved by the smallest adjusted p-value and then the smallest
genomic start coordinate.

`start` and `end` must be zero-based, half-open BED coordinates describing one
base. `browser_start_1based` must equal `start + 1`.

The strand is read from `PASid`, for example:

```text
chr1:+:123456
chr2:-:987654
```

The coordinate inside `PASid` is not used.

### GENCODE GTF

The GTF must contain transcript features with:

```text
tag "Ensembl_canonical"
```

The example APAlyzer table contains gene symbols, so the script normally
matches `gene_symbol` to GTF `gene_name`, chromosome, and strand. If the APA
table contains a nonblank `gene_id`, that stable identifier is preferred.

The first run creates a reusable gffutils SQLite database in the output
directory. Later runs reuse it, which saves substantial time on an HPC.

Use `--rebuild-gtf-db` after changing the GTF. Alternatively, place the
database at a specific shared path with `--gtf-db`.

### Sample sheet

The sample sheet is tab-separated:

```text
sample_id	condition	replicate	role	pair_id	bigwig
control_1	control	1	control	1	/path/control_1.bw
control_2	control	2	control	2	/path/control_2.bw
drugA_1	drugA	1	treatment	1	/path/drugA_1.bw
drugA_2	drugA	2	treatment	2	/path/drugA_2.bw
drugB_1	drugB	1	treatment	1	/path/drugB_1.bw
drugB_2	drugB	2	treatment	2	/path/drugB_2.bw
```

`role` must be `control` or `treatment`. Multiple treatment conditions are
allowed.

Relative BigWig paths are resolved relative to the sample sheet's directory.

#### Paired controls

When `pair_id` is filled, every treatment is compared with the one control
having the same `pair_id`.

#### Unpaired controls

For an unpaired run, leave `pair_id` blank for every sample. Each treatment
replicate is compared with the mean of all control BigWigs.

The script prints a warning because the plotted SD then reflects treatment
replicate variation but does not fully preserve control-replicate variation.

Do not mix paired and unpaired samples in one run.

### BigWig requirements

BigWigs must contain:

- nonnegative signal
- linear-scale signal
- comparable normalization across samples

Signed tracks or tracks that are already log2 transformed cannot be used with
the ratio calculation in this script.

Missing BigWig positions are treated as zero.

## Pseudocount

The script calculates this value for every gene and normalized position:

```text
log2((treatment signal + pseudocount) /
     (control signal + pseudocount))
```

`--pseudocount` is required because an appropriate value depends on the scale
of the normalized BigWigs. For example, `0.1` may be reasonable for one signal
scale and inappropriate for another. Record the chosen value in analysis
notes and consider checking that conclusions are not driven by near-zero
signal.

## Running the Script

```bash
python dipa_metagene.py \
    --apa-results apalyzer_results.tsv \
    --gtf gencode.annotation.gtf \
    --samples samples.tsv \
    --output-dir metagene_results \
    --pseudocount 0.1 \
    --adjusted-pvalue 0.05 \
    --bins-per-side 100 \
    --min-side-bp 100 \
    --threads 8
```

`--min-side-bp` must be at least as large as `--bins-per-side`. This guarantees
that every normalized bin contains at least one real nucleotide. No
interpolation or duplicated bases are used.

Run `python dipa_metagene.py --help` for all options.

## Threading on an HPC

BigWig extraction is parallelized by sample with `ThreadPoolExecutor`.
Each thread opens its own BigWig handle; handles are never shared.

The maximum useful thread count is the number of BigWig samples. Requesting 32
threads for 6 samples will use 6 threads.

Set `--threads 0` to detect CPU allocations from these variables:

- `SLURM_CPUS_PER_TASK`
- `PBS_NP`
- `NSLOTS`

If none is defined, the script uses the CPUs visible to Python, capped at the
number of samples.

### SLURM example

```bash
#!/bin/bash
#SBATCH --job-name=dipa_metagene
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=dipa_metagene.%j.log

module load python/3.11
source /path/to/dipa-metagene-env/bin/activate

python /path/to/dipa_metagene.py \
    --apa-results /path/to/apalyzer_results.tsv \
    --gtf /path/to/gencode.annotation.gtf \
    --samples /path/to/samples.tsv \
    --output-dir /path/to/metagene_results \
    --pseudocount 0.1 \
    --threads "$SLURM_CPUS_PER_TASK"
```

If the BigWigs are on a heavily shared filesystem, more threads may eventually
stop improving speed because storage bandwidth becomes the limit.

## Aggregation

The script gives every retained gene equal weight:

1. Extract and bin signal separately for each gene and sample.
2. Calculate gene-level treatment-versus-control log2FC.
3. Average genes within each treatment replicate.
4. Average replicate metagene profiles for the treatment line.
5. Calculate sample standard deviation across replicate profiles at every bin.

With one replicate, the SD is missing and no shaded band is drawn.

The same final gene set is used for every treatment and sample.

## Outputs

The output directory contains:

### Plots

- `combined_metagene.png`: 300 DPI image
- `combined_metagene.pdf`: vector plot

### Tables

- `selected_up_genes.tsv`: final genes, canonical transcripts, dIPA context,
  removed exon information, and cDNA lengths
- `excluded_genes.tsv`: every filtered or excluded row and the reason
- `per_gene_log2fc.tsv`: gene-level log2FC for each treatment replicate and bin
- `replicate_profiles.tsv`: gene-averaged profile for every treatment replicate
- `treatment_summary.tsv`: treatment mean and per-bin replicate SD

### Reproducibility

- `run_parameters.json`: paths, parameters, package versions, thread count,
  treatment counts, and filtering totals
- `*.gffutils.db`: reusable GTF database unless `--gtf-db` points elsewhere

## Common Errors

### No matching Ensembl canonical transcript

Check that:

- the GTF is GENCODE and includes `Ensembl_canonical` tags
- gene symbols correspond to GTF `gene_name`
- chromosome names agree
- APA and GTF genome assemblies agree

### dIPA outside the canonical transcript

The APA coordinates and GTF may use different genome assemblies, or the dIPA
may have been assigned through a different transcript model.

### Insufficient cDNA

After splitting at an intron or removing the dIPA-containing exon, one side has
fewer than `--min-side-bp` retained exonic bases. The gene is excluded rather
than stretched to artificial resolution.

### Negative BigWig value

The BigWig is signed or transformed. Use comparable, nonnegative,
linear-normalized signal tracks.

### Chromosome missing from BigWig

Check naming such as `chr1` versus `1` and confirm genome assemblies.

### Existing GTF database is wrong

Rerun with:

```bash
--rebuild-gtf-db
```

## Tests

After installing dependencies:

```bash
pytest -q
```

The tests create small temporary GTF, APA, sample-sheet, and BigWig files.
They do not require real experimental data.
