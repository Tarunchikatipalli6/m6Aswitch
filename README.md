# m6Aswitch

An R package for integrating m6A (N6-methyladenosine) predictions with isoform switching events.

## Overview

m6Aswitch joins two datasets:

- **m6A sites** from m6Anet — transcript, position, modification probability
- **Isoform switches** from IsoformSwitchAnalyzeR — gene, isoform pair, FDR, dIF

and reports what happens to m6A sites when a gene switches isoforms.

Because m6A sites are compared at genomic coordinates rather than transcript positions, results remain valid when the two isoforms differ in length or exon composition.

## Installation

```r
devtools::install_github("Tarunchikatipalli6/m6Aswitch")
```

## Quick start

```r
library(m6Aswitch)

GTF <- "gencode.v49.annotation.gtf"

# Parse m6A sites from each condition
m6a_wt  <- parse_m6anet("wt.csv",  probability_threshold = 0.9)
m6a_mut <- parse_m6anet("mut.csv", probability_threshold = 0.9)

# Tag conditions and combine
m6a_wt[,  condition := "WT"]
m6a_mut[, condition := "MUT"]
combined <- rbind(m6a_wt, m6a_mut)

# Parse isoform switches
iso <- parse_isoform_switch("switch_pairs.txt", fdr_threshold = 0.05)

# Lift to genomic coordinates and annotate
gr  <- lift_m6a_to_genomic(combined, GTF)
res <- annotate_m6a_switches_genomic(gr, iso)
res <- annotate_drach(res, m6a_wt, m6a_mut)
res <- classify_lost_mechanism(res, GTF)

# Results
table(res$isoform_status)
table(res$m6a_fate)

export_annotated_switches(res, "results")
plot_m6aswitch_results(res, iso, "plots", "WT vs MUT")
```

## Two classifications

The package reports two independent columns, answering different questions.

**`isoform_status`** — which isoform of the pair carries the site

| Value | Meaning |
|---|---|
| `ISOFORM_A_ONLY` | site is on isoform A only |
| `ISOFORM_B_ONLY` | site is on isoform B only |
| `IN_BOTH_ISOFORMS` | the genomic position exists on both |

This reflects transcript structure. Always computed.

**`m6a_fate`** — which condition detected the site

| Value | Meaning |
|---|---|
| `LOST` | detected in condition 1 only |
| `GAINED` | detected in condition 2 only |
| `RETAINED` | detected in both |

This reflects methylation state. Computed only when a `condition` column is supplied to `lift_m6a_to_genomic()`.

A site can be present in both conditions while sitting on only one isoform — `ISOFORM_B_ONLY` with `m6a_fate = RETAINED`. Read both columns.

## Functions

### `parse_m6anet()`

```r
parse_m6anet(m6anet_file,
             probability_threshold = 0.9,
             transcript_col = "transcript_id",
             position_col   = "position",
             prob_col       = "probability",
             keep_extra     = NULL)
```

Reads m6Anet output and standardises column names. Converts m6Anet's 0-based positions to 1-based.

For m6Anet's own output naming:

```r
parse_m6anet("data.site_proba.csv",
             transcript_col = "ensembl_transcript_id",
             position_col   = "transcript_position",
             prob_col       = "probability_modified")
```

`keep_extra` retains additional columns such as `n_reads`.

**Returns:** `transcript_id`, `position`, `probability`, `kmer`

### `parse_isoform_switch()`

```r
parse_isoform_switch(iso_switch_file,
                     fdr_threshold = 0.05,
                     gene_col   = "geneID",
                     iso_a_col  = "isoformID_A",
                     iso_b_col  = "isoformID_B",
                     fdr_col    = "isoform_switch_q_value")
```

Reads switch pairs, filters by FDR, and deduplicates keeping the most significant pair per gene-isoform triple.

Prints which condition each isoform is favoured in, since that direction depends on factor ordering upstream.

**Returns:** `gene_id`, `isoform_a`, `isoform_b`, `fdr`, `condition_1`, `condition_2`, `dif`

### `lift_m6a_to_genomic()`

```r
lift_m6a_to_genomic(m6a_sites, gtf_file)
```

Converts transcript coordinates to genomic coordinates using exon structure from the GTF.

This step is required because transcript position 150 in two isoforms with different exon composition maps to different genomic locations. Comparing transcript positions directly produces incorrect calls.

Additional columns, including `condition`, are carried through.

**Returns:** GRanges with `transcript_id`, `transcript_position`, `probability`

### `annotate_m6a_switches_genomic()`

```r
annotate_m6a_switches_genomic(m6a_sites_gr, iso_switches)
```

Compares m6A sites at each genomic locus across the isoform pair, and across conditions when condition information is available.

**Returns:** data.table with `isoform_status`, `m6a_fate`, `probability_a`, `probability_b`, `probability_c1`, `probability_c2`, plus gene, isoform, coordinate and switch statistics.

### `annotate_drach()`

```r
annotate_drach(m6a_switches, m6a_condition_a, m6a_condition_b)
```

Adds a `drach_motif` column using the kmer reported by m6Anet, so no transcript sequence file is needed.

DRACH is **D** (A/G/U) · **R** (A/G) · **A** · **C** · **H** (A/C/U) — 18 possible 5-mers.

m6Anet only evaluates DRACH positions, so this column will be uniformly TRUE for m6Anet input. It becomes informative with tools that test other contexts.

### `classify_lost_mechanism()`

```r
classify_lost_mechanism(m6a_switches, gtf_file)
```

For sites present on only one isoform, determines whether the other isoform contains that genomic position.

| Value | Type | Meaning |
|---|---|---|
| `A_ONLY_EXON_SKIPPED` | structural | position absent from isoform B |
| `A_ONLY_UNMETHYLATED` | regulatory | position present in B, not methylated |
| `B_ONLY_EXON_INCLUDED` | structural | position absent from isoform A |
| `B_ONLY_NEW_METHYLATION` | regulatory | position present in A, methylated only in B |

Structural differences are explained by splicing. Regulatory differences are not.

**Adds:** `isoform_mechanism`, `target_isoform_has_exon`

### `export_annotated_switches()`

```r
export_annotated_switches(m6a_switches,
                          output_prefix = "m6aswitch_results",
                          format        = "both",
                          add_igv_track = TRUE,
                          color_by      = "isoform_status")
```

Writes CSV and/or BED9 with IGV-compatible colouring.

### `plot_m6aswitch_results()`

```r
plot_m6aswitch_results(classified, switch_pairs, output_dir,
                       analysis_name = "Analysis",
                       color_by      = "isoform_status")
```

Generates five PDFs: classification distribution, mechanism breakdown, top-genes heatmap, dIF versus Δprobability volcano, and per-gene track plots.

`plot_m6a_switches()` returns a single plot for interactive use.

## Input formats

**m6Anet output**

```
transcript_id,transcript_position,n_reads,probability_modified,kmer,mod_ratio
ENST00000054950.4,1481,42,0.8121,GGACT,0.31
```

Positions are 0-based; the package converts them.

**Switch pairs**

```
geneID	isoformID_A	isoformID_B	condition_1	condition_2	isoform_switch_q_value	dIF
GENE1	ENST00000001.1	ENST00000002.1	WT	MUT	0.001	0.35
```

Built from `isoformSwitchTestDEXSeq()` output by pairing each gene's up- and down-regulated isoform:

```r
iso_dt <- as.data.table(iso_dexseq$isoformFeatures)
sig <- iso_dt[!is.na(isoform_switch_q_value) &
              isoform_switch_q_value < 0.05 & abs(dIF) >= 0.1]

pairs <- list()
for (g in unique(sig$gene_id)) {
  sub  <- sig[gene_id == g]
  up   <- sub[dIF > 0]
  down <- sub[dIF < 0]
  if (nrow(up) > 0 && nrow(down) > 0) {
    pairs[[length(pairs) + 1]] <- data.table(
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
fwrite(rbindlist(pairs), "switch_pairs.txt", sep = "\t")
```

**GTF** — must match the annotation release used for alignment.

## Notes on interpretation

**Isoform A and B.** IsoformSwitchAnalyzeR defines `dIF = IF(condition_2) − IF(condition_1)`. Isoform A is the one favoured in condition 1, isoform B in condition 2. Condition order follows factor level, which may not match expectation. `parse_isoform_switch()` prints the direction; set it explicitly with `factor(condition, levels = c("WT", "MUT"))` in the design matrix if needed.

**Detection versus absence.** m6Anet reports a position only when read coverage is sufficient. A position missing from one condition may be unmethylated, or may never have been evaluated. When library depths differ substantially, restrict the analysis to positions reported in both conditions:

```r
a_all <- parse_m6anet("cond_a.csv", probability_threshold = 0)
b_all <- parse_m6anet("cond_b.csv", probability_threshold = 0)

testable <- merge(unique(a_all[, .(transcript_id, position)]),
                  unique(b_all[, .(transcript_id, position)]),
                  by = c("transcript_id", "position"))

m6a_a <- merge(parse_m6anet("cond_a.csv", 0.9), testable,
               by = c("transcript_id", "position"))
m6a_b <- merge(parse_m6anet("cond_b.csv", 0.9), testable,
               by = c("transcript_id", "position"))
```

**Isoform length.** A longer isoform contains more DRACH positions and will tend to carry more detected sites. Check transcript lengths before interpreting an asymmetry in `isoform_status`.

## Threshold selection

| Threshold | Use |
|---|---|
| 0.9 | Primary results |
| 0.8 | Moderate |
| 0.5 | Sensitivity analysis |
| 0.0 | All evaluated positions |

Running both a strict and a permissive threshold and reporting them separately is recommended.

## Testing

```r
devtools::test()
```

## Dependencies

`data.table`, `GenomicRanges`, `GenomicFeatures`, `txdbmaker`, `IRanges`, `S4Vectors`, `Biostrings`, `methods`, `ggplot2`, `tidyr`, `stringr`, `readr`, `igraph`

## Citation

```
Chikatipalli, T. (2026). m6Aswitch: Integration of m6A modifications with
isoform switching. GitHub: Tarunchikatipalli6/m6Aswitch
```

Methodology:

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

## Author

Tarun Chikatipalli — [GitHub](https://github.com/Tarunchikatipalli6)

## License

MIT — see LICENSE.

## Issues

[github.com/Tarunchikatipalli6/m6Aswitch/issues](https://github.com/Tarunchikatipalli6/m6Aswitch/issues)
