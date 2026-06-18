# Comprehensive Demo Script for m6Aswitch Package
# Demonstrates dual-threshold analysis: HIGH-CONFIDENCE (0.9) + SENSITIVITY (0.5)
# 
# This workflow aligns with publication best practices:
# - PRIMARY RESULTS: High-confidence threshold (0.9)
# - SUPPLEMENTARY: Sensitivity analysis (0.5)

# ============================================================================
# SETUP
# ============================================================================

library(m6Aswitch)
library(data.table)

# Set working directory to the sample_data folder
sample_data_dir <- "~/Desktop/m6aswitch/sample_data"

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("m6Aswitch: Dual-Threshold Analysis Workflow\n")
cat("════════════════════════════════════════════════════════════════\n")

# ============================================================================
# STEP 1: PARSE INPUT DATA (Once, used for both analyses)
# ============================================================================

cat("\n========== STEP 1: Parse Input Data ==========\n")

# Parse isoform switches (same for both thresholds)
iso_switch_file <- file.path(sample_data_dir, "sample_isoform_switches.txt")
iso_switches <- parse_isoform_switch(
  iso_switch_file = iso_switch_file,
  fdr_threshold = 0.05,
  gene_col = "geneID",
  iso_a_col = "isoformID_A",
  iso_b_col = "isoformID_B",
  fdr_col = "isoform_switch_q_value"
)

cat("✓ Loaded", nrow(iso_switches), "isoform switches\n")

# Load isoform sequences (same for both thresholds)
seq_file <- file.path(sample_data_dir, "sample_isoform_sequences.csv")
iso_sequences <- fread(seq_file)

cat("✓ Loaded sequences for", nrow(iso_sequences), "isoforms\n")

# ============================================================================
# ANALYSIS 1: HIGH-CONFIDENCE m6A (probability ≥ 0.9) - PRIMARY RESULTS
# ============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS 1: HIGH-CONFIDENCE m6A DETECTION (threshold ≥ 0.9)\n")
cat("Purpose: Primary publication-quality results\n")
cat("════════════════════════════════════════════════════════════════\n")

# STEP 1A: Parse m6A with HIGH threshold
cat("\n--- Step 1A: Parse m6A predictions (probability ≥ 0.9) ---\n")

m6a_file <- file.path(sample_data_dir, "sample_m6a_predictions.csv")
m6a_high <- parse_m6anet(
  m6anet_file = m6a_file,
  probability_threshold = 0.9,
  transcript_col = "transcript_id",
  position_col = "position",
  prob_col = "probability"
)

cat("✓ Identified", nrow(m6a_high), "high-confidence m6A sites\n")
cat("  (probability ≥ 0.9)\n")
print(head(m6a_high, 3))

# STEP 1B: Annotate m6A changes (HIGH confidence)
cat("\n--- Step 1B: Annotate m6A changes (HIGH confidence) ---\n")

m6a_switches_high <- annotate_m6a_switches(
  m6a_sites = m6a_high,
  iso_switches = iso_switches,
  iso_sequences = iso_sequences
)

cat("✓ Identified", nrow(m6a_switches_high), "m6A-isoform interactions\n")

# STEP 1C: Add DRACH motif annotation (HIGH confidence)
cat("\n--- Step 1C: Add DRACH motif validation (HIGH confidence) ---\n")

m6a_switches_high_final <- annotate_drach(
  m6a_switches = m6a_switches_high,
  sequences = iso_sequences
)

cat("✓ Added DRACH motif annotations\n")

# STEP 1D: Summary statistics (HIGH confidence)
cat("\n--- HIGH-CONFIDENCE RESULTS SUMMARY ---\n")
cat("Total m6A sites (prob ≥ 0.9):", nrow(m6a_high), "\n")
cat("Total m6A-switch interactions:", nrow(m6a_switches_high_final), "\n")

if (nrow(m6a_switches_high_final) > 0) {
  cat("\nm6A Fate Distribution (HIGH confidence):\n")
  fate_dist_high <- table(m6a_switches_high_final$m6a_fate)
  print(fate_dist_high)

  cat("\nGenes affected (HIGH confidence):\n")
  genes_high <- unique(m6a_switches_high_final$gene_id)
  print(genes_high)
  cat("Number of unique genes:", length(genes_high), "\n")

  cat("\nDRACH motif validation (HIGH confidence):\n")
  drach_dist_high <- table(m6a_switches_high_final$drach_motif)
  print(drach_dist_high)
}

# ============================================================================
# ANALYSIS 2: SENSITIVITY m6A (probability ≥ 0.5) - SUPPLEMENTARY RESULTS
# ============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS 2: SENSITIVITY m6A DETECTION (threshold ≥ 0.5)\n")
cat("Purpose: Exploratory analysis, supplementary materials\n")
cat("════════════════════════════════════════════════════════════════\n")

# STEP 2A: Parse m6A with LENIENT threshold
cat("\n--- Step 2A: Parse m6A predictions (probability ≥ 0.5) ---\n")

m6a_broad <- parse_m6anet(
  m6anet_file = m6a_file,
  probability_threshold = 0.5,
  transcript_col = "transcript_id",
  position_col = "position",
  prob_col = "probability"
)

cat("✓ Identified", nrow(m6a_broad), "sensitivity m6A sites\n")
cat("  (probability ≥ 0.5)\n")
print(head(m6a_broad, 3))

# STEP 2B: Annotate m6A changes (SENSITIVITY)
cat("\n--- Step 2B: Annotate m6A changes (SENSITIVITY) ---\n")

m6a_switches_broad <- annotate_m6a_switches(
  m6a_sites = m6a_broad,
  iso_switches = iso_switches,
  iso_sequences = iso_sequences
)

cat("✓ Identified", nrow(m6a_switches_broad), "m6A-isoform interactions\n")

# STEP 2C: Add DRACH motif annotation (SENSITIVITY)
cat("\n--- Step 2C: Add DRACH motif validation (SENSITIVITY) ---\n")

m6a_switches_broad_final <- annotate_drach(
  m6a_switches = m6a_switches_broad,
  sequences = iso_sequences
)

cat("✓ Added DRACH motif annotations\n")

# STEP 2D: Summary statistics (SENSITIVITY)
cat("\n--- SENSITIVITY RESULTS SUMMARY ---\n")
cat("Total m6A sites (prob ≥ 0.5):", nrow(m6a_broad), "\n")
cat("Total m6A-switch interactions:", nrow(m6a_switches_broad_final), "\n")

if (nrow(m6a_switches_broad_final) > 0) {
  cat("\nm6A Fate Distribution (SENSITIVITY):\n")
  fate_dist_broad <- table(m6a_switches_broad_final$m6a_fate)
  print(fate_dist_broad)

  cat("\nGenes affected (SENSITIVITY):\n")
  genes_broad <- unique(m6a_switches_broad_final$gene_id)
  print(genes_broad)
  cat("Number of unique genes:", length(genes_broad), "\n")

  cat("\nDRACH motif validation (SENSITIVITY):\n")
  drach_dist_broad <- table(m6a_switches_broad_final$drach_motif)
  print(drach_dist_broad)
}

# ============================================================================
# COMPARATIVE ANALYSIS: HIGH vs SENSITIVITY
# ============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("COMPARATIVE ANALYSIS: HIGH-CONFIDENCE (0.9) vs SENSITIVITY (0.5)\n")
cat("════════════════════════════════════════════════════════════════\n")

n_high <- nrow(m6a_switches_high_final)
n_broad <- nrow(m6a_switches_broad_final)
n_additional <- n_broad - n_high
pct_additional <- if (n_broad > 0) round(100 * n_additional / n_broad, 1) else 0

cat("\nSummary Comparison:\n")
cat("  High-confidence (≥0.9):  ", n_high, "m6A-switch interactions\n")
cat("  Sensitivity (≥0.5):       ", n_broad, "m6A-switch interactions\n")
cat("  Additional (0.5-0.9):     ", n_additional, 
    paste0("(", pct_additional, "% of total)\n"))

genes_high <- unique(m6a_switches_high_final$gene_id)
genes_broad <- unique(m6a_switches_broad_final$gene_id)
genes_shared <- length(intersect(genes_high, genes_broad))
genes_only_broad <- setdiff(genes_broad, genes_high)

cat("\nGene-level Comparison:\n")
cat("  Genes in high-confidence:     ", length(genes_high), "\n")
cat("  Genes in sensitivity:         ", length(genes_broad), "\n")
cat("  Genes in both:                ", genes_shared, "\n")
cat("  Genes only in sensitivity:    ", length(genes_only_broad), "\n")

# ============================================================================
# EXPORT RESULTS
# ============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("EXPORT RESULTS\n")
cat("════════════════════════════════════════════════════════════════\n")

# Export HIGH confidence results
cat("\n--- Exporting HIGH-confidence results ---\n")
output_high_csv <- file.path(sample_data_dir, "m6a_switch_results_HIGH_CONFIDENCE.csv")
fwrite(m6a_switches_high_final, file = output_high_csv)
cat("✓ High-confidence: ", output_high_csv, "\n")

# Export SENSITIVITY results
cat("\n--- Exporting SENSITIVITY results ---\n")
output_broad_csv <- file.path(sample_data_dir, "m6a_switch_results_SENSITIVITY.csv")
fwrite(m6a_switches_broad_final, file = output_broad_csv)
cat("✓ Sensitivity:      ", output_broad_csv, "\n")

# Export comparison summary
cat("\n--- Exporting comparison summary ---\n")
output_summary <- file.path(sample_data_dir, "m6a_analysis_summary.txt")
con <- file(output_summary, "w")
writeLines(c(
  "m6Aswitch: Dual-Threshold Analysis Summary",
  "==========================================",
  "",
  "HIGH-CONFIDENCE ANALYSIS (probability ≥ 0.9)",
  paste("  m6A sites detected:          ", nrow(m6a_high)),
  paste("  m6A-switch interactions:     ", nrow(m6a_switches_high_final)),
  paste("  Unique genes affected:       ", length(genes_high)),
  "",
  "SENSITIVITY ANALYSIS (probability ≥ 0.5)",
  paste("  m6A sites detected:          ", nrow(m6a_broad)),
  paste("  m6A-switch interactions:     ", nrow(m6a_switches_broad_final)),
  paste("  Unique genes affected:       ", length(genes_broad)),
  "",
  "COMPARATIVE METRICS",
  paste("  Additional interactions (0.5-0.9): ", n_additional),
  paste("  Additional genes:                    ", length(genes_only_broad)),
  "",
  "OUTPUT FILES",
  paste("  High-confidence: ", output_high_csv),
  paste("  Sensitivity:     ", output_broad_csv)
), con)
close(con)
cat("✓ Summary:         ", output_summary, "\n")

# ============================================================================
# FINAL MESSAGE
# ============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS COMPLETE\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("\nFor your publication:\n")
cat("  PRIMARY RESULTS:    Use HIGH-CONFIDENCE (threshold ≥ 0.9)\n")
cat("  SUPPLEMENTARY:      Use SENSITIVITY (threshold ≥ 0.5)\n")
cat("  MANUSCRIPT TEXT:    'We identified X isoform switches with\n")
cat("                       high-confidence m6A changes (probability ≥ 0.9,\n")
cat("                       n=", n_high, "). A sensitivity analysis at\n")
cat("                       probability ≥ 0.5 identified Y additional\n")
cat("                       interactions (n=", n_additional, ",\n")
cat("                       Supplementary Table X).'\n")
cat("\n")
