# Python dIPA Metagene Plot: Implementation Brief

## Purpose

Write a beginner-friendly Python script that creates a combined metagene plot
centered conceptually on differential intronic polyadenylation sites (dIPAs).

The genes and dIPAs come from an APAlyzer results table. Only significant
events in the `UP` direction are used. The signal comes from normalized
per-replicate BigWig files. Each treatment is plotted as log2 fold-change
relative to its control, with a shaded standard deviation across biological
replicates at each normalized position.

The signal must represent spliced cDNA/exonic sequence, not the full genomic
gene body.

## Number-One Coding Requirement

Readability for a beginner is more important than concision.

Use a single, straightforward Python script with clearly labeled sections,
descriptive variable names, and short comments that explain biological or
coordinate-system decisions. Avoid clever shortcuts, dense comprehensions,
advanced object-oriented design, classes, decorators, and unnecessary
abstractions.

Only create functions for operations that are genuinely repeated, such as
extracting and orienting BigWig values for many exons or binning many cDNA
profiles. Do not turn every small step into a separate function. The main
workflow should read from top to bottom in a clear, sequential order.

Prefer explicit intermediate variables over nested expressions. Error messages
must explain what is wrong and how the user can fix it.

## Deliverables

Create:

1. `dipa_metagene.py`
2. `requirements.txt`
3. `README.md`
4. `tests/test_dipa_metagene.py`

The script should be runnable from the command line.

Example:

```bash
python dipa_metagene.py \
    --apa-results apalyzer_results.tsv \
    --gtf gencode.annotation.gtf \
    --samples samples.tsv \
    --output-dir metagene_results \
    --pseudocount 0.1 \
    --adjusted-pvalue 0.05 \
    --bins-per-side 100 \
    --min-side-bp 100
```

Do not assume that `0.1` is universally appropriate. The README must explain
that the pseudocount should be selected based on the scale of the normalized
BigWigs.

## Python Dependencies

Use familiar, well-supported packages:

- `pandas` for tabular input and output
- `numpy` for signal arrays and calculations
- `gffutils` for structured GTF parsing
- `pyBigWig` for BigWig access
- `matplotlib` for plotting
- `seaborn` only for a colorblind-friendly palette, if useful
- `pytest` for tests

Avoid adding dependencies unless they clearly simplify an important task.

## Input 1: APAlyzer Results

The APAlyzer table is tab-separated and contains columns similar to:

```text
gene_symbol
PASid
RED
pvalue
p_adj
APAreg
APAreg_original
significance_for_plot
chr
start
end
browser_start_1based
browser_end_1based
browser_region
pas_coordinate_match
```

Validate that all required columns are present. Stop with a readable error that
lists any missing columns.

### Meaning of RED

`RED` is APAlyzer's Relative Expression Difference. A positive value indicates
increased relative IPA usage in the treatment comparison, while a negative
value indicates decreased usage.

Use `RED` only for:

1. Confirming that an `UP` event has a positive direction.
2. Choosing one representative dIPA when a gene has multiple significant UP
   sites.

Do not use `RED` as the plotted y-axis value. The plotted values must be
calculated from the BigWigs.

### Filter Significant UP Events

Keep rows satisfying all of these conditions:

```text
p_adj <= adjusted-pvalue threshold
APAreg == "UP"
RED > 0
pas_coordinate_match == TRUE
```

Also check `APAreg_original`. If it does not equal `UP`, exclude the row and
record the reason as `APAreg_and_APAreg_original_disagree`.

Do not use `significance_for_plot` to decide whether a row is significant. It
is a plotting helper column, not the authoritative filtering field.

Handle common text encodings of true values, such as `TRUE`, `True`, and
`true`, without making the input case-sensitive.

### Select One dIPA Per Gene

The metagene should give each gene equal weight. Therefore, use one
representative dIPA per significant UP gene.

When a gene has multiple qualifying rows:

1. Select the row with the largest positive `RED`.
2. If RED values tie, select the row with the smallest `p_adj`.
3. If both tie, select the row with the smallest genomic `start` coordinate so
   the choice is deterministic.

Retain the unselected rows in the QC output with reason
`additional_significant_PAS_for_gene`.

In the supplied example, filtering gives 28 UP PAS rows representing 16 genes
before transcript and cDNA-length validation.

## Coordinate Conventions

Coordinate handling must be explicit because off-by-one errors would move the
dIPA or assign it to the wrong exon.

- `start` and `end` are BED-style coordinates for BigWig access:
  zero-based, half-open `[start, end)`.
- `browser_start_1based` is the one-based dIPA position used when comparing
  with GTF coordinates.
- GTF feature coordinates are one-based and inclusive.
- Use `chr`, `start`, and `end` as the authoritative mapped genomic location.
- Do not use the coordinate embedded inside `PASid`; it differs from the mapped
  coordinates in the provided example.
- Parse the strand (`+` or `-`) from `PASid`.

Validate:

```text
browser_start_1based == start + 1
end == start + 1
```

If either validation fails, exclude the row and explain the mismatch in the QC
table.

## Input 2: GENCODE GTF

Use `gffutils` to parse the GTF. Do not parse the entire attributes column with
fragile string splitting.

### Match the Gene

The example APAlyzer table contains `gene_symbol` rather than a stable Ensembl
gene ID. Match `gene_symbol` to the GTF `gene_name` attribute.

If a future input table contains a usable `gene_id` column, prefer `gene_id`
over `gene_symbol`.

Use chromosome and strand to resolve or validate the match. If the symbol is
missing or ambiguous after chromosome and strand validation, exclude the gene
and record a specific reason.

### Select the Canonical Transcript

For each selected gene, use only the transcript whose GTF attributes include:

```text
tag "Ensembl_canonical"
```

The canonical flag selects the transcript model used to define the exons.

Do not silently fall back to the longest transcript or another isoform.

Exclude and report genes when:

- no `Ensembl_canonical` transcript is present
- more than one canonical transcript remains after gene, chromosome, and
  strand matching
- the canonical transcript has fewer than two exons
- the dIPA lies outside the canonical transcript span

## Determine the dIPA Context

Sort the canonical transcript's exons in transcriptional 5-prime to 3-prime
order.

- For a plus-strand transcript, this is increasing genomic coordinate.
- For a minus-strand transcript, this is decreasing genomic coordinate.

Determine whether the one-based dIPA coordinate falls inside a canonical exon
or between two adjacent canonical exons.

### Intronic dIPA

If the dIPA lies between two adjacent canonical exons:

- The upstream cDNA side contains all exons through the exon immediately
  before the dIPA-containing intron in transcript order.
- The downstream cDNA side contains all exons beginning with the exon
  immediately after that intron.
- Do not include any part of the intron.
- Record the context as `intronic`.

### Exonic dIPA

If the dIPA lies inside a canonical exon:

- Remove the entire containing exon from the cDNA profile.
- The upstream cDNA side contains all exons before the removed exon in
  transcript order.
- The downstream cDNA side contains all exons after the removed exon.
- Do not query or plot BigWig signal from the removed exon.
- Record the context as `exonic_exon_removed`.
- Record the removed exon ID and its genomic coordinates in the selected-gene
  output.

If the dIPA is exactly on an exon coordinate, treat it as exonic because GTF
exon coordinates are inclusive.

## Minimum cDNA Sequence on Each Side

The default plot uses 100 normalized bins before and 100 normalized bins after
the conceptual dIPA center.

These are normalized bins, not necessarily 100 genomic bases. However, do not
stretch fewer than 100 real exonic bases into 100 bins because that would imply
artificial resolution.

Calculate the total included exonic length separately for the upstream and
downstream sides.

Require:

```text
upstream_exonic_bp >= min-side-bp
downstream_exonic_bp >= min-side-bp
```

The default `--min-side-bp` is 100.

Exclude genes that fail with one of these reasons:

```text
insufficient_upstream_cdna
insufficient_downstream_cdna
insufficient_cdna_on_both_sides
```

This rule also excludes a dIPA in the first or last exon when removal of that
exon leaves one side empty.

## Input 3: Sample Sheet

Use a tab-separated sample sheet with these columns:

```text
sample_id
condition
replicate
role
pair_id
bigwig
```

Example:

```text
sample_id	condition	replicate	role	pair_id	bigwig
control_1	control	1	control	1	/path/control_1.bw
control_2	control	2	control	2	/path/control_2.bw
drugA_1	drugA	1	treatment	1	/path/drugA_1.bw
drugA_2	drugA	2	treatment	2	/path/drugA_2.bw
drugB_1	drugB	1	treatment	1	/path/drugB_1.bw
drugB_2	drugB	2	treatment	2	/path/drugB_2.bw
```

Requirements:

- `sample_id` values must be unique.
- `role` must be either `control` or `treatment`.
- At least one control and one treatment must be present.
- Multiple treatment conditions are allowed.
- Every BigWig path must exist and be readable.
- BigWig chromosome names must be compatible with the GTF and APAlyzer table.
- BigWig values must be nonnegative and on a comparable linear normalized
  scale. Log-transformed or signed BigWigs are not valid for the ratio formula
  below.

### Control Pairing

When `pair_id` is present:

- Match each treatment replicate to the control with the same `pair_id`.
- The same paired control may be used for different treatment conditions.
- Stop with a readable error if a treatment has no matching control or more
  than one matching control.

If pairing information is absent or blank for all samples:

- Calculate the mean control signal for each gene and bin.
- Compare each treatment replicate with that mean control profile.
- Print a warning that the SD then reflects variation among treatment
  replicates but does not fully preserve control-replicate variation.

Do not mix paired and unpaired treatment samples in the same run.

## Extract BigWig Signal

Use `pyBigWig.values(chromosome, start, end, numpy=True)` to obtain base-level
signal for each retained exon.

For each exon:

1. Query with zero-based, half-open coordinates.
2. Replace missing `NaN` values with zero by default.
3. Keep values in genomic order for plus-strand transcripts.
4. Reverse the values for minus-strand transcripts.

After orienting individual exons, concatenate them in transcript 5-prime to
3-prime order.

Process one gene and sample at a time rather than loading every chromosome into
memory.

If a BigWig lacks a required chromosome, exclude the gene from the entire run
so every condition uses the same gene set. Record the affected chromosome and
sample in QC.

If any extracted value is negative, stop with a clear explanation that
log2-ratio normalization requires nonnegative linear signal.

## Normalize Each cDNA Side into Bins

Use the same number of bins for every gene so profiles can be averaged.

Default:

```text
100 upstream bins
100 downstream bins
200 total plotted bins
```

Because each retained side has at least as many bases as bins, use
`numpy.array_split()` to divide the concatenated upstream values into equal or
nearly equal groups. Calculate the arithmetic mean of each group. Repeat for
the downstream side.

Do not interpolate signal and do not duplicate bases to fill bins.

Concatenate the 100 upstream means and 100 downstream means into one
200-position profile.

## Treatment-versus-Control Normalization

The user must provide a positive pseudocount through `--pseudocount`.

For every retained gene, normalized position, treatment replicate, and matched
control, calculate:

```python
log2_fold_change = numpy.log2(
    (treatment_signal + pseudocount)
    / (control_signal + pseudocount)
)
```

If controls are unpaired, replace `control_signal` with the mean control signal
for that gene and normalized position.

Perform the log2FC calculation at the gene-and-bin level before averaging
genes. Do not first average all genes in treatment and control and then take a
ratio.

## Aggregate Genes and Replicates

Use the same final gene set for all samples and treatments.

Each gene contributes one representative dIPA and receives equal weight.

For each treatment replicate:

1. Calculate the log2FC profile for every retained gene.
2. Average the gene-level log2FC values at each normalized position.
3. Save this as the replicate metagene profile.

For each treatment condition:

1. Calculate the mean of its replicate metagene profiles at every position.
2. Calculate the sample standard deviation across replicate profiles using
   `ddof=1`.
3. If only one replicate exists, do not draw an SD band and store SD as
   missing rather than zero.

The shaded region is therefore per-normalized-position variability across
biological replicate metagene profiles. Call it "per-bin SD" in code and
documentation, since the x-axis is normalized rather than literal genomic
bases.

## Plot

Make one combined plot with:

- one line per treatment condition
- a colorblind-friendly color for each treatment
- a matching transparent band showing mean plus or minus one per-bin SD
- a horizontal line at log2FC zero
- a vertical line at the conceptual dIPA center
- x-axis labels `0%`, `dIPA`, and `100%`
- y-axis label `log2 fold-change vs control`
- title `Metagene signal around UP-regulated dIPAs`
- a legend listing treatment conditions
- the final number of genes and replicate counts in a caption or subtitle

The upstream and downstream sides should occupy equal visual width even though
their original cDNA lengths differ among genes.

Use a clear white background and readable font sizes. Save:

```text
combined_metagene.png
combined_metagene.pdf
```

Save the PNG at 300 DPI.

## Tabular Outputs

Write these files into `--output-dir`.

### `selected_up_genes.tsv`

Include at least:

```text
gene_symbol
gene_id
transcript_id
chromosome
strand
RED
p_adj
dipa_start_0based
dipa_position_1based
dipa_context
removed_exon_id
removed_exon_start
removed_exon_end
upstream_exonic_bp
downstream_exonic_bp
```

Use blank values for removed-exon fields when the dIPA is intronic.

### `excluded_genes.tsv`

Include the original APAlyzer identifiers plus:

```text
exclusion_stage
exclusion_reason
details
```

Do not silently discard rows.

### `per_gene_log2fc.tsv`

Store gene-level profiles in long format:

```text
gene_symbol
condition
sample_id
pair_id
position_index
region
normalized_position
log2fc
```

`region` should be `upstream` or `downstream`.

### `replicate_profiles.tsv`

```text
condition
sample_id
pair_id
position_index
region
normalized_position
mean_gene_log2fc
number_of_genes
```

### `treatment_summary.tsv`

```text
condition
position_index
region
normalized_position
mean_log2fc
sd_log2fc
number_of_replicates
number_of_genes
```

### `run_parameters.json`

Record:

- all command-line arguments
- package versions
- input paths
- filtering counts
- final gene count
- treatment and replicate counts
- whether controls were paired
- pseudocount
- date and time

## Console Summary

At the end, print a short, readable summary such as:

```text
APAlyzer rows read: 29
Significant UP PAS rows: 28
Unique significant UP genes: 16
Genes with a canonical transcript: 15
Intronic dIPAs retained: 10
Exonic dIPAs retained after exon removal: 3
Genes excluded for short cDNA sides: 2
Final genes used in every profile: 13
Treatments plotted: drugA, drugB
Results written to: /path/metagene_results
```

Use the real counts from the run.

## Error Handling

Fail early with beginner-friendly messages for:

- missing input files
- missing required columns
- malformed numeric fields
- invalid adjusted p-value threshold
- zero or negative pseudocount
- zero bins or minimum side length
- duplicate sample IDs
- missing or ambiguous controls
- mismatched chromosome naming
- unreadable or invalid BigWigs
- negative BigWig values
- no significant UP events
- no genes remaining after annotation or length filtering

An error should identify the relevant gene, sample, path, or column whenever
possible.

## Testing Requirements

Use small synthetic inputs so tests are fast and understandable.

At minimum, test:

1. Significant UP filtering.
2. Exclusion of `DN`, nonsignificant, mismatched-direction, and failed
   coordinate rows.
3. Selection of the largest RED per gene and deterministic tie breaking.
4. BED-to-GTF coordinate conversion.
5. Selection of the `Ensembl_canonical` transcript.
6. Plus-strand intronic dIPA splitting.
7. Minus-strand intronic dIPA splitting and signal reversal.
8. Exonic dIPA removal on both strands.
9. Exclusion when removal of the first or last exon leaves an empty side.
10. Exclusion at 99 bp and retention at exactly 100 bp with default settings.
11. Binning without interpolation or empty bins.
12. Correct paired log2FC calculation.
13. Unpaired mean-control fallback.
14. Equal gene weighting.
15. Mean and sample SD across replicate profiles.
16. Multiple treatments sharing the same control set.
17. Rejection of negative BigWig values.

When practical, create tiny synthetic BigWigs with `pyBigWig` inside temporary
test directories.

## README Requirements

Write the README for someone who is comfortable running command-line tools but
is not an experienced Python programmer.

Explain:

- what the script calculates
- what RED means and how it is used
- why the plotted y-axis is independently calculated from BigWigs
- expected coordinate systems
- how `Ensembl_canonical` is used
- what happens to an exon containing a dIPA
- why at least 100 exonic bases are required on each side by default
- the sample-sheet format
- paired versus unpaired controls
- how to choose the pseudocount
- all output files
- common errors and how to fix them

Include a complete example command and a small example sample sheet.

## Acceptance Criteria

The implementation is complete when:

1. A user can run one documented Python command using an APAlyzer table, a
   GENCODE GTF, a sample sheet, and BigWigs.
2. Only significant UP genes are included.
3. One highest-RED dIPA is used per gene.
4. Only the `Ensembl_canonical` transcript defines the cDNA profile.
5. Intronic dIPAs split the profile at the host intron.
6. Exonic dIPAs cause the entire containing exon to be omitted.
7. Both cDNA sides contain at least the configured minimum number of exonic
   bases.
8. Minus-strand profiles are correctly oriented 5-prime to 3-prime.
9. Every treatment is shown as BigWig-derived log2FC relative to control.
10. Each treatment line is the replicate mean and its band is per-bin replicate
    SD.
11. Every exclusion is visible in QC output.
12. The script, comments, README, and tests prioritize clarity for a beginner.
