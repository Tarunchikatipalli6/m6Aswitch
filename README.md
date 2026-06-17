# m6Aswitch

An R package for analyzing dynamic m6A (N6-methyladenosine) modifications across isoform switches in RNA.

## Overview

**m6Aswitch** integrates m6A methylation predictions with isoform switching events to identify functional changes in RNA modifications during alternative splicing. The package helps researchers:

- Parse m6A predictions from m6AnetAnalyzer
- Process isoform switching data from IsoformSwitchAnalyzeR
- Detect m6A changes (LOST, GAINED, RETAINED) across isoform pairs
- Annotate m6A sites with DRACH motif information
- Characterize functional implications of m6A switches

## Installation

```r
# Install from GitHub
devtools::install_github("Tarunchikatipalli6/m6Aswitch")
```

## Quick Start

```r
library(m6Aswitch)
library(data.table)

# Parse m6A predictions
m6a_sites <- parse_m6anet("m6anet_predictions.csv", probability_threshold = 0.5)

# Parse isoform switches
iso_switches <- parse_isoform_switch("isoform_switches.txt", fdr_threshold = 0.05)

# Load isoform sequences
iso_sequences <- fread("isoform_sequences.csv")

# Annotate m6A changes across switches
results <- annotate_m6a_switches(m6a_sites, iso_switches, iso_sequences)

# Add DRACH motif annotations
results <- annotate_drach(results, iso_sequences)

# View results
View(results)
```

## Functions

### Parsing Functions

#### `parse_m6anet()`
Reads m6A site predictions from m6AnetAnalyzer output and returns standardized data.

```r
m6a_sites <- parse_m6anet(
  m6anet_file = "predictions.csv",
  probability_threshold = 0.5,
  transcript_col = "transcript_id",
  position_col = "position",
  prob_col = "probability"
)
```

**Output Columns:**
- `transcript_id`: Transcript identifier
- `position`: Genomic position of m6A site
- `probability`: Prediction probability (0-1)
- `gene_id`: Gene identifier (optional)
- `strand`: Strand orientation (optional)

---

#### `parse_isoform_switch()`
Reads isoform switching results from IsoformSwitchAnalyzeR.

```r
iso_switches <- parse_isoform_switch(
  iso_switch_file = "switches.txt",
  fdr_threshold = 0.05,
  gene_col = "geneID",
  iso_a_col = "isoformID_A",
  iso_b_col = "isoformID_B",
  fdr_col = "isoform_switch_q_value"
)
```

**Output Columns:**
- `gene_id`: Gene identifier
- `isoform_a`: First isoform ID
- `isoform_b`: Second isoform ID
- `fdr`: FDR-adjusted p-value
- Additional columns from input file

---

### Annotation Functions

#### `annotate_m6a_switches()`
Integrates m6A sites with isoform switches to detect functional changes.

```r
m6a_switches <- annotate_m6a_switches(
  m6a_sites,
  iso_switches,
  iso_sequences
)
```

**Output Columns:**
- All columns from `iso_switches`
- `position`: m6A site position
- `m6a_in_isoform_a`: Logical, is m6A present in isoform A?
- `m6a_in_isoform_b`: Logical, is m6A present in isoform B?
- `m6a_fate`: Character - "LOST", "GAINED", or "RETAINED"
- `probability_a`: m6A prediction probability in isoform A
- `probability_b`: m6A prediction probability in isoform B

**m6A Fate Categories:**
- **LOST**: m6A present in isoform A but not in isoform B
- **GAINED**: m6A present in isoform B but not in isoform A
- **RETAINED**: m6A present in both isoforms

---

#### `annotate_drach()`
Adds DRACH motif annotations to m6A sites.

```r
results <- annotate_drach(
  m6a_switches,
  sequences = iso_sequences
)
```

**DRACH Pattern:**
- **D** = A/G/U (purine or uracil)
- **R** = A/G (purine)
- **A** = A (adenosine - m6A target)
- **C** = C (cytosine)
- **H** = A/C/U (any except G)

**Output Columns:**
- `drach_motif`: Logical - does the site match DRACH?
- `drach_motif_a`: Logical - DRACH status in isoform A
- `drach_motif_b`: Logical - DRACH status in isoform B

---

#### `find_motif()`
Detects DRACH motif at a specific position in a sequence.

```r
has_drach <- find_motif(
  sequence = "ACGAACATCG",
  position = 5,
  context_bp = 2
)
```

---

### Utility Functions

#### `export_annotated_switches()`
Export results to CSV format.

```r
export_annotated_switches(
  m6a_switches = results,
  filename = "m6a_switch_results.csv"
)
```

---

## Sample Data

The package includes sample data files in the `sample_data/` directory:

- `sample_m6a_predictions.csv` - m6A site predictions
- `sample_isoform_switches.txt` - Isoform switching events
- `sample_isoform_sequences.csv` - RNA sequences

### Run the Demo

```r
source("~/Desktop/m6aswitch/demo_workflow.R")
```

This will process the sample data and generate:
- Summary statistics
- m6A fate distribution
- Affected genes
- Export results to CSV

---

## Testing

The package includes comprehensive unit tests (33 tests, all passing).

```r
# Run all tests
devtools::test(pkg = "~/Desktop/m6aswitch")

# Run specific test file
devtools::test_file("~/Desktop/m6aswitch/tests/testthat/test-find_motif.R")
```

See `TESTING.md` for detailed testing guide.

---

## Input File Formats

### m6A Predictions (CSV/TSV)
```
transcript_id,position,probability,gene_id,strand
ENST00000001,100,0.95,ENSG00001,+
ENST00000002,150,0.88,ENSG00001,+
```

### Isoform Switches (TSV)
```
geneID	isoformID_A	isoformID_B	isoform_switch_q_value	dIF
ENSG00001	ENST00000001	ENST00000002	0.001	-0.35
```

### Isoform Sequences (CSV)
```
isoform_id,sequence
ENST00000001,ACGAACATCGGAACACCGGGAACAAACGAAC...
ENST00000002,ACGAAGATCGGAACACCGGGAACAAACGAAC...
```

---

## Coordinate System and Genomic Lifting

### The Problem

m6AnetAnalyzer reports m6A sites in **transcript-level coordinates**. This means
position 245 in `ENST00000001` and position 245 in `ENST00000002` refer to
completely different genomic loci. Comparing raw transcript positions across
isoforms produces spurious LOST/GAINED/RETAINED calls.

### The Solution: GTF-Aware Coordinate Lifting

Use `lift_m6a_to_genomic()` to convert transcript coordinates to genomic
coordinates before annotating switches. Then use
`annotate_m6a_switches_genomic()` to compare sites at the **same genomic
location** across isoforms.

```r
library(m6Aswitch)
library(data.table)

# 1. Parse inputs (as usual)
m6a_sites     <- parse_m6anet("m6anet_predictions.csv")
iso_switches  <- parse_isoform_switch("isoform_switches.txt")
iso_sequences <- fread("isoform_sequences.csv")

# 2. Lift transcript coordinates -> genomic coordinates (requires GTF)
genomic_sites <- lift_m6a_to_genomic(m6a_sites, gtf_file = "genome.gtf")
# Returns a GRanges object; genomic_sites$transcript_id and
# genomic_sites$probability carry over from m6a_sites

# 3. Annotate switches using genomic coordinates
m6a_annotated <- annotate_m6a_switches_genomic(
  genomic_sites, iso_switches, iso_sequences
)

# 4. Optional: add DRACH motif annotation (uses transcript sequences)
m6a_annotated_drach <- annotate_drach(m6a_annotated, iso_sequences)

# 5. Visualize
plot_m6a_switches(m6a_annotated, plot_type = "summary")
```

### Before / After Comparison

| Approach | Comparison basis | Scientific validity |
|----------|-----------------|---------------------|
| `annotate_m6a_switches()` | Transcript position number | ❌ Position 245 means different things per isoform |
| `annotate_m6a_switches_genomic()` | Genomic coordinates (chr:start-end) | ✅ Same locus compared across isoforms |

---

## Example Workflow Output

```
========== SUMMARY STATISTICS ==========

Total m6A sites analyzed: 14 
Total isoform switches: 5 
Total m6A-switch interactions: 21 

m6A Fate Distribution:
  GAINED     LOST RETAINED 
       9       10        2 

Genes affected: ENSG00001, ENSG00002, ENSG00003
```

---

## Citation

If you use m6Aswitch in your research, please cite:

```
Chikatipalli, T. (2026). m6Aswitch: Analysis of dynamic m6A modifications 
across isoform switches. GitHub: Tarunchikatipalli6/m6Aswitch
```

---

## Dependencies

- `data.table` - Fast data manipulation
- `tidyr` - Data reshaping
- `Biostrings` - Biological string handling
- `GenomicRanges` - Genomic data structures
- `GenomicFeatures` - GTF/TxDb handling and coordinate mapping
- `IRanges` - Range objects
- `S4Vectors` - S4 vector infrastructure
- `methods` - R methods infrastructure
- `igraph` - Network analysis
- `ggplot2` - Visualization
- `stringr` - String operations
- `readr` - File reading

---

## Features

✅ Parse m6A predictions from m6AnetAnalyzer  
✅ Process isoform switches from IsoformSwitchAnalyzeR  
✅ Detect m6A changes (LOST/GAINED/RETAINED)  
✅ GTF-aware coordinate lifting (transcript → genomic)  
✅ Scientifically valid cross-isoform comparison via genomic coordinates  
✅ DRACH motif annotation  
✅ Comprehensive unit tests  
✅ Sample data and demo workflow  
✅ Easy-to-use API  

---

## Documentation

- `TESTING.md` - Complete testing guide
- `demo_workflow.R` - Example analysis
- `sample_data/` - Sample input files

---

## Author

Tarun Chikatipalli  
[GitHub](https://github.com/Tarunchikatipalli6)

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

---

## Contact

For questions, issues, or contributions, please open an issue on [GitHub](https://github.com/Tarunchikatipalli6/m6Aswitch/issues).
