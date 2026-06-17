# Demo Script for m6Aswitch Package
# This script demonstrates the complete workflow using sample data

# ============================================================================
# SETUP
# ============================================================================

library(m6Aswitch)
library(data.table)

# Set working directory to the sample_data folder
sample_data_dir <- "~/Desktop/m6aswitch/sample_data"

# ============================================================================
# STEP 1: PARSE m6A PREDICTIONS
# ============================================================================

cat("\n========== STEP 1: Parse m6A Predictions ==========\n")

m6a_file <- file.path(sample_data_dir, "sample_m6a_predictions.csv")
m6a_sites <- parse_m6anet(
  m6anet_file = m6a_file,
  probability_threshold = 0.80,
  transcript_col = "transcript_id",
  position_col = "position",
  prob_col = "probability"
)

cat("✓ Loaded", nrow(m6a_sites), "m6A sites\n")
print("First few m6A sites:")
print(head(m6a_sites, 3))

# ============================================================================
# STEP 2: PARSE ISOFORM SWITCHES
# ============================================================================

cat("\n========== STEP 2: Parse Isoform Switches ==========\n")

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
print("Isoform switches:")
print(iso_switches)

# ============================================================================
# STEP 3: LOAD ISOFORM SEQUENCES
# ============================================================================

cat("\n========== STEP 3: Load Isoform Sequences ==========\n")

seq_file <- file.path(sample_data_dir, "sample_isoform_sequences.csv")
iso_sequences <- fread(seq_file)

cat("✓ Loaded sequences for", nrow(iso_sequences), "isoforms\n")
print("Sequence data (first 3):")
print(head(iso_sequences, 3))

# ============================================================================
# STEP 4: ANNOTATE m6A CHANGES ACROSS ISOFORM SWITCHES
# ============================================================================

cat("\n========== STEP 4: Annotate m6A Changes ==========\n")

m6a_switches <- annotate_m6a_switches(
  m6a_sites = m6a_sites,
  iso_switches = iso_switches,
  iso_sequences = iso_sequences
)

cat("✓ Annotated", nrow(m6a_switches), "m6A-isoform switch interactions\n")
print("Annotated results (first 5 rows):")
print(head(m6a_switches, 5))

# ============================================================================
# STEP 5: ADD DRACH MOTIF ANNOTATIONS
# ============================================================================

cat("\n========== STEP 5: Annotate DRACH Motifs ==========\n")

m6a_switches_final <- annotate_drach(
  m6a_switches = m6a_switches,
  sequences = iso_sequences
)

cat("✓ Added DRACH motif annotations\n")
print("Final annotated results (first 5 rows):")
print(head(m6a_switches_final, 5))

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n========== SUMMARY STATISTICS ==========\n")

cat("\nTotal m6A sites analyzed:", nrow(m6a_sites), "\n")
cat("Total isoform switches:", nrow(iso_switches), "\n")
cat("Total m6A-switch interactions:", nrow(m6a_switches_final), "\n")

if (nrow(m6a_switches_final) > 0) {
  cat("\nm6A Fate Distribution:\n")
  fate_dist <- table(m6a_switches_final$m6a_fate)
  print(fate_dist)
  
  cat("\nGenes affected:\n")
  genes <- unique(m6a_switches_final$gene_id)
  print(genes)
  
  cat("\nNumber of unique genes:", length(genes), "\n")
}

# ============================================================================
# EXPORT RESULTS (Optional)
# ============================================================================

cat("\n========== EXPORT RESULTS ==========\n")

output_file <- file.path(sample_data_dir, "m6a_switch_results.csv")
fwrite(m6a_switches_final, file = output_file)
cat("✓ Results exported to:", output_file, "\n")

# ============================================================================
# VISUALIZATION (Optional - if you want to add plots later)
# ============================================================================

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("All files have been processed successfully!\n")
cat("Check the results in: sample_data/m6a_switch_results.csv\n")
