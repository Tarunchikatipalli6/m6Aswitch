# m6Aswitch

An R package for integrating m6A (N6-methyladenosine) predictions with isoform switching events.

---

## What the package answers

You have two datasets:

- **m6A sites** from m6Anet — transcript, position, probability
- **Isoform switches** from IsoformSwitchAnalyzeR — gene, isoform A, isoform B, FDR, dIF

m6Aswitch joins them and answers **two distinct questions**:

| Question | Column | Type |
|---|---|---|
| Do the two isoforms of a switch carry different m6A sites? | `isoform_status` | structural |
| Does methylation at a position differ between conditions? | `m6a_fate` | regulatory |

These are frequently conflated. A site can be detected in **both conditions** while sitting on **only one isoform** of the pair — in which case `isoform_status` is `ISOFORM_B_ONLY` but `m6a_fate` is `RETAINED`.

Read both columns. They are not interchangeable.

---

## Installation

```r
devtools::install_github("Tarunchikatipalli6/m6Aswitch")
```

---

## Quick start

```r
library(m6Aswitch)

GTF <- "gencode.v49.annotation.gtf"

# ── 1. Find positions m6Anet evaluated in BOTH conditions ───────────────────
a_all <- parse_m6anet("cond_a.csv", probability_threshold = 0,
                      transcript_col = "ensembl_transcript_id",
                      position_col   = "transcript_position",
                      prob_col       = "probability_modified")
b_all <- parse_m6anet("cond_b.csv", probability_threshold = 0, ...)

testable <- merge(unique(a_all[, .(transcript_id, position)]),
                  unique(b_all[, .(transcript_id, position)]),
                  by = c("transcript_id", "position"))

# ── 2. Parse at threshold, restrict to testable positions ───────────────────
m6a_a <- merge(parse_m6anet("cond_a.csv", 0.9, ...), testable,
               by = c("transcript_id", "position"))
m6a_b <- merge(parse_m6anet("cond_b.csv", 0.9, ...), testable,
               by = c("transcript_id", "position"))

# ── 3. Tag conditions and combine — do NOT deduplicate ──────────────────────
m6a_a[, condition := "WT"]
m6a_b[, condition := "MUT"]
combined <- rbind(m6a_a, m6a_b)

# ── 4. Parse switches — check the printed A/B direction ─────────────────────
iso <- parse_isoform_switch("switch_pairs.txt", fdr_threshold = 0.05)

# ── 5. Lift, annotate, classify ─────────────────────────────────────────────
gr  <- lift_m6a_to_genomic(combined, GTF)
res <- annotate_m6a_switches_genomic(gr, iso)
res <- annotate_drach(res, m6a_a, m6a_b)
res <- classify_lost_mechanism(res, GTF)

# ── 6. Results ──────────────────────────────────────────────────────────────
table(res$isoform_status)                 # structural
table(res$m6a_fate, useNA = "ifany")      # regulatory
table(res$isoform_status, res$m6a_fate)   # are they independent?

export_annotated_switches(res, "results")
plot_m6aswitch_results(res, iso, "plots", "My Analysis")
```

---

## Two things that will bite you

### 1. Absence is not the same as untested

m6Anet reports a position **only when coverage is sufficient**. A position missing from one condition's output may be unmethylated — or simply never evaluated.

Treating these as equivalent produces spurious gains when one library is shallower. In one real analysis:

```
Astrocyte (1 replicate) :  5,095 sites
Tumour (3 pooled)       : 28,644 sites     5.6x

Result: 80.6% "GAINED"

Diagnostic: of 4,612 GAINED sites,
            4,611 were NEVER MEASURED in astrocyte
```

After restricting to positions evaluated in both:

```
Astrocyte : 4,300     Tumour : 4,141
Condition-level: RETAINED 574, LOST 136, GAINED 108
```

**Always build a testability filter** — steps 1 and 2 in the Quick Start.

### 2. Which isoform is A?

IsoformSwitchAnalyzeR defines `dIF = IF(condition_2) − IF(condition_1)`, and switch pairs are normally built as:

```r
isoformID_A = down$isoform_id[1]   # dIF < 0, favoured in condition_1
isoformID_B = up$isoform_id[1]     # dIF > 0, favoured in condition_2
```

Condition order comes from factor level, which may not match your intuition — `"IDH_MUT"` sorts before `"IDH_WT"`, making the mutant condition_1.

`parse_isoform_switch()` prints the direction:

```
Parsed switch_pairs.txt: 53 switch pair(s) at FDR <= 0.05
  Comparison: IDH_R132H (condition_1) vs IDH_WT (condition_2)
  isoform_a is favoured in IDH_R132H
  isoform_b is favoured in IDH_WT
```

To control it, set the level explicitly before running IsoformSwitchAnalyzeR:

```r
design$condition <- factor(design$condition, levels = c("WT", "MUT"))
```

---

## Function reference

### `parse_m6anet()`

```r
parse_m6anet(m6anet_file,
             probability_threshold = 0.9,
             transcript_col = "transcript_id",
             position_col   = "position",
             prob_col       = "probability",
             keep_extra     = NULL)
```

Reads m6Anet output. Converts **0-based** m6Anet positions to **1-based**.

The conversion is verified against transcript sequence — at the raw position the base is usually G; at position + 1 it is always A, and the reported kmer matches only after the shift.

For m6Anet's own file naming:

```r
parse_m6anet("data.site_proba.csv",
             transcript_col = "ensembl_transcript_id",
             position_col   = "transcript_position",
             prob_col       = "probability_modified")
```

Use `keep_extra` to carry through additional columns (`n_reads`, `condition`).

| Threshold | Use |
|---|---|
| 0.9 | Primary, publication |
| 0.8 | Moderate |
| 0.5 | Sensitivity analysis |
| 0.0 | Every evaluated position — for the testability filter |

**Returns:** `transcript_id`, `position`, `probability`, `kmer`

---

### `parse_isoform_switch()`

```r
parse_isoform_switch(iso_switch_file, fdr_threshold = 0.05,
                     gene_col = "geneID", iso_a_col = "isoformID_A",
                     iso_b_col = "isoformID_B",
                     fdr_col = "isoform_switch_q_value")
```

Reads switch pairs, filters by FDR, deduplicates keeping the most significant, and prints the A/B direction.

**Returns:** `gene_id`, `isoform_a`, `isoform_b`, `fdr`, `condition_1`, `condition_2`, `dif`

---

### `lift_m6a_to_genomic()`

```r
lift_m6a_to_genomic(m6a_sites, gtf_file)
```

Converts transcript coordinates to genomic coordinates using exon structure from the GTF.

**Why this is necessary:** transcript position 150 in two isoforms with different exon composition maps to different genomic locations.

```
Isoform A exons:
  Exon1  chr10:1,000-1,100    (transcript positions 1-100)
  Exon2  chr10:5,000-5,200    (transcript positions 101-300)

Transcript position 150 = 50 bases into Exon2
  = chr10:5,049                    CORRECT

Using transcript boundaries alone:
  = 1,000 + 150 = chr10:1,150      WRONG (inside the intron)
```

Extra columns — including `condition` — are forwarded to the output GRanges.

**Returns:** GRanges with `transcript_id`, `transcript_position`, `probability`, plus forwarded columns.

---

### `annotate_m6a_switches_genomic()`

```r
annotate_m6a_switches_genomic(m6a_sites_gr, iso_switches)
```

Produces two independent classifications.

**`isoform_status` — structural**

| Value | Meaning |
|---|---|
| `ISOFORM_A_ONLY` | site is on isoform A only |
| `ISOFORM_B_ONLY` | site is on isoform B only |
| `IN_BOTH_ISOFORMS` | the genomic position exists on both |

Always computed. Reflects transcript structure — a longer isoform carries more sites.

**`m6a_fate` — regulatory**

| Value | Meaning |
|---|---|
| `LOST` | detected in condition_1 only |
| `GAINED` | detected in condition_2 only |
| `RETAINED` | detected in both |

`NA` unless a `condition` column was supplied to `lift_m6a_to_genomic()`. A warning is issued when it is missing, and again if the condition labels don't match those in `iso_switches`.

---

### `annotate_drach()`

```r
annotate_drach(m6a_switches, m6a_condition_a, m6a_condition_b)
```

Adds `drach_motif` using the kmer m6Anet reports — no transcript sequence file needed.

DRACH: **D** = A/G/U · **R** = A/G · **A** = A · **C** = C · **H** = A/C/U → 18 possible 5-mers.

**Note:** m6Anet only evaluates DRACH positions, so all its calls are DRACH by construction. For m6Anet input this column confirms rather than filters. It becomes informative with tools that test non-DRACH positions.

---

### `classify_lost_mechanism()`

```r
classify_lost_mechanism(m6a_switches, gtf_file)
```

For sites present on only one isoform, asks whether the *other* isoform contains that genomic position.

| Value | Type | Meaning |
|---|---|---|
| `A_ONLY_EXON_SKIPPED` | structural | position absent from isoform B |
| `A_ONLY_UNMETHYLATED` | regulatory | position present in B, not methylated |
| `B_ONLY_EXON_INCLUDED` | structural | position absent from isoform A |
| `B_ONLY_NEW_METHYLATION` | regulatory | position present in A, methylated only in B |

The **regulatory** classes cannot be explained by splicing.

Adds `isoform_mechanism` and `target_isoform_has_exon`.

---

### `export_annotated_switches()`

```r
export_annotated_switches(m6a_switches, output_prefix = "results",
                          format = "both", add_igv_track = TRUE,
                          color_by = "isoform_status")
```

CSV and/or BED9 with IGV colouring. `color_by` selects which classification the BED encodes.

---

### `plot_m6aswitch_results()`

```r
plot_m6aswitch_results(classified, switch_pairs, output_dir,
                       analysis_name = "Analysis",
                       color_by = "isoform_status")
```

Five PDFs: distribution, mechanism breakdown, top-genes heatmap, dIF-vs-Δprobability volcano, and per-gene track plots.

`color_by = "m6a_fate"` colours by condition instead — errors clearly if condition information wasn't supplied.

---

## Interpreting results

A realistic output from 35 switch pairs carrying m6A:

```
Isoform status:
  ISOFORM_B_ONLY      78
  ISOFORM_A_ONLY      36
  IN_BOTH_ISOFORMS    32

Mechanism:
  Structural (splicing):
    A_ONLY_EXON_SKIPPED       2
    B_ONLY_EXON_INCLUDED      0
  Regulatory (position in both isoforms):
    A_ONLY_UNMETHYLATED      34
    B_ONLY_NEW_METHYLATION   78
```

**Reading it:** 112 of 114 differences (98%) occur at positions **both isoforms contain**. Only 2 are explained by exon skipping. So isoform switching is not relocating m6A sites — the differences reflect modification state.

**What not to read:** "78 sites gained methylation in condition B." `isoform_status` describes which transcript carries the site. Use `m6a_fate` for the condition comparison.

**Why B often dominates:** isoform B may simply be longer or better covered — more sequence, more DRACH positions, more detected sites. Check transcript lengths before interpreting an asymmetry.

---

## Input formats

**m6Anet** (`data.site_proba.csv`)

```
transcript_id,transcript_position,n_reads,probability_modified,kmer,mod_ratio
ENST00000054950.4,1481,42,0.8121,GGACT,0.31
```

Positions are 0-based; the package converts them.

**Switch pairs** (TSV)

```
geneID	isoformID_A	isoformID_B	condition_1	condition_2	isoform_switch_q_value	dIF
GENE1	ENST00000001.1	ENST00000002.1	WT	MUT	0.001	0.35
```

Built from `isoformSwitchTestDEXSeq()` output:

```r
iso_dt <- as.data.table(iso_dexseq$isoformFeatures)
sig <- iso_dt[!is.na(isoform_switch_q_value) &
              isoform_switch_q_value < 0.05 & abs(dIF) >= 0.1]

pairs_list <- list()
for (g in unique(sig$gene_id)) {
  sub  <- sig[gene_id == g]
  up   <- sub[dIF > 0]
  down <- sub[dIF < 0]
  if (nrow(up) > 0 && nrow(down) > 0) {
    pairs_list[[length(pairs_list) + 1]] <- data.table(
      geneID      = g,
      isoformID_A = down$isoform_id[1],
      isoformID_B = up$isoform_id[1],
      condition_1 = down$condition_1[1],
      condition_2 = down$condition_2[1],
      isoform_switch_q_value = min(sub$isoform_switch_q_value),
      dIF = up$dIF[1]
    )
  }
}
fwrite(rbindlist(pairs_list), "switch_pairs.txt", sep = "\t")
```

**GTF** — must be the same annotation release the reads were aligned to. Mixing versions puts coordinates in the wrong place.

---

## Validation

The coordinate handling has been verified end to end:

| Check | Result |
|---|---|
| Round-trip (transcript → genome → transcript) | 178/178 exact |
| Manual exon-walk verification | matches |
| Plus-strand correlation | median +0.999 |
| Minus-strand correlation | median −0.992 |
| Sites landing within exons | 178/178 |
| Base at converted position | 100/100 are A |
| Reported kmer matches sequence | 100/100 |
| Fate vs presence flags | 0 violations |
| Mechanism vs independent recompute | 0 disagreements |

---

## Dual-threshold analysis

Running a strict and a permissive threshold, reported separately, is good practice.

```r
# Primary
m6a_high <- parse_m6anet("predictions.csv", probability_threshold = 0.9)

# Sensitivity
m6a_broad <- parse_m6anet("predictions.csv", probability_threshold = 0.5)
```

| | 0.9 | 0.5 |
|---|---|---|
| Confidence | high | moderate |
| False positives | low | moderate |
| Sensitivity | lower | higher |
| Use | main results | supplementary |

Both should be applied **after** the testability filter, not instead of it. The two address different problems: the threshold controls call confidence; the filter controls whether a position was evaluated at all.

---

## Testing

```r
devtools::test()
```

---

## Dependencies

`data.table` · `GenomicRanges` · `GenomicFeatures` · `txdbmaker` · `IRanges` · `S4Vectors` · `Biostrings` · `methods` · `ggplot2` · `tidyr` · `stringr` · `readr` · `igraph`

---

## Citation

```
Chikatipalli, T. (2026). m6Aswitch: Integration of m6A modifications with
isoform switching. GitHub: Tarunchikatipalli6/m6Aswitch
```

Methodology references:

```
Vitting-Seerup, K. & Sandelin, A. (2019). IsoformSwitchAnalyzeR: analysis of
changes in genome-wide patterns of alternative splicing and its functional
consequences. Bioinformatics 35(21):4469-4471.

Hendra, C., et al. (2022). Detection of m6A from direct RNA sequencing using
a multiple instance learning framework. Nature Methods 19, 1590-1598.

Liu, H., Begik, O., Lucas, M. C., et al. (2023). Accurate detection of m6A
RNA modifications in the eukaryotic transcriptome with SCARLET.
Nature Biotechnology 41, 896-905.
```

---

## Author

Tarun Chikatipalli · [GitHub](https://github.com/Tarunchikatipalli6)

## License

MIT — see LICENSE.

## Issues

[github.com/Tarunchikatipalli6/m6Aswitch/issues](https://github.com/Tarunchikatipalli6/m6Aswitch/issues)
