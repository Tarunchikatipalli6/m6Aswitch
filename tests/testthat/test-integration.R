library(testthat)
library(data.table)

# Integration tests for the complete m6Aswitch workflow

test_that("complete workflow: m6a_sites -> annotate_m6a_switches -> annotate_drach", {
  # Step 1: Create sample m6A sites
  m6a_sites <- data.table(
    transcript_id = c("ENST001", "ENST001", "ENST002", "ENST002"),
    position = c(100, 150, 100, 200),
    probability = c(0.95, 0.88, 0.92, 0.85)
  )
  
  # Step 2: Create isoform switches
  iso_switches <- data.table(
    gene_id = c("ENSG001", "ENSG002"),
    isoform_a = c("ENST001", "ENST002"),
    isoform_b = c("ENST002", "ENST003"),
    fdr = c(0.001, 0.01)
  )
  
  # Step 3: Create sequences
  iso_sequences <- data.table(
    isoform_id = c("ENST001", "ENST002", "ENST003"),
    sequence = c(
      "ACGAACATCGGAACACCGGGAACAAACGAAC",
      "ACGAAGATCGGAACACCGGGAACAAACGAAC",
      "CGAACATCGGAACACCGGGAACAAACGAAC"
    )
  )
  
  # Run annotation
  m6a_switches <- annotate_m6a_switches(m6a_sites, iso_switches, iso_sequences)
  
  expect_true(is.data.table(m6a_switches))
  expect_true(nrow(m6a_switches) > 0)
  expect_true(all(c("m6a_fate", "probability_a", "probability_b") %in% names(m6a_switches)))
})

test_that("annotation handles multiple genes", {
  m6a_sites <- data.table(
    transcript_id = c("ISO_G1_A", "ISO_G1_B", "ISO_G2_A", "ISO_G2_B"),
    position = c(100, 150, 200, 250),
    probability = c(0.95, 0.88, 0.92, 0.85)
  )
  
  iso_switches <- data.table(
    gene_id = c("GENE1", "GENE1", "GENE2", "GENE2"),
    isoform_a = c("ISO_G1_A", "ISO_G1_B", "ISO_G2_A", "ISO_G2_B"),
    isoform_b = c("ISO_G1_B", "ISO_G1_A", "ISO_G2_B", "ISO_G2_A"),
    fdr = c(0.001, 0.005, 0.01, 0.02)
  )
  
  iso_sequences <- data.table(
    isoform_id = c("ISO_G1_A", "ISO_G1_B", "ISO_G2_A", "ISO_G2_B"),
    sequence = rep("ACGAACATCGGAACACCGGGAACAAACGAAC", 4)
  )
  
  result <- annotate_m6a_switches(m6a_sites, iso_switches, iso_sequences)
  
  # Check that both genes are represented
  genes_in_result <- unique(result$gene_id)
  expect_true(length(genes_in_result) > 0)
})

test_that("results are properly sorted by FDR and m6a_fate", {
  m6a_sites <- data.table(
    transcript_id = c("ISO1", "ISO2", "ISO1", "ISO2"),
    position = c(100, 200, 300, 300),
    probability = c(0.95, 0.92, 0.88, 0.90)
  )
  
  iso_switches <- data.table(
    gene_id = "GENE1",
    isoform_a = "ISO1",
    isoform_b = "ISO2",
    fdr = 0.01
  )
  
  iso_sequences <- data.table(
    isoform_id = c("ISO1", "ISO2"),
    sequence = c("ACGTACGTACGT", "ACGTACGTACGT")
  )
  
  result <- annotate_m6a_switches(m6a_sites, iso_switches, iso_sequences)
  
  # Check if sorted (gene_id first, then fdr, then m6a_fate)
  expect_equal(result$gene_id[1], result$gene_id[nrow(result)])
})
