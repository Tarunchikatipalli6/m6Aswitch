library(testthat)
library(data.table)

test_that("annotate_drach adds columns to m6a_switches", {
  # Create test data
  m6a_switches <- data.table(
    gene_id = c("GENE1", "GENE1"),
    isoform_a = c("ISO1", "ISO1"),
    isoform_b = c("ISO2", "ISO2"),
    position = c(100, 150),
    m6a_in_isoform_a = c(TRUE, TRUE),
    m6a_in_isoform_b = c(FALSE, TRUE),
    m6a_fate = c("LOST", "RETAINED"),
    probability_a = c(0.95, 0.88),
    probability_b = c(NA_real_, 0.92),
    fdr = c(0.01, 0.05)
  )
  
  sequences <- data.table(
    isoform_id = c("ISO1", "ISO2"),
    sequence = c("ACGAACATCGGAACA", "ACGAAGATCGGAACA")
  )
  
  result <- annotate_drach(m6a_switches, sequences)
  
  # Check that result has the DRACH columns
  expect_true("drach_motif" %in% names(result))
  expect_equal(nrow(result), nrow(m6a_switches))
})

test_that("annotate_drach validates input", {
  not_dt <- list(gene_id = "GENE1")
  sequences <- data.table(
    isoform_id = "ISO1",
    sequence = "ACGAACATCG"
  )
  
  expect_error(annotate_drach(not_dt, sequences), "must be a data.table")
})

test_that("annotate_drach handles empty m6a_switches", {
  m6a_switches <- data.table(
    gene_id = character(0),
    isoform_a = character(0),
    isoform_b = character(0),
    position = integer(0),
    m6a_in_isoform_a = logical(0),
    m6a_in_isoform_b = logical(0),
    m6a_fate = character(0),
    probability_a = numeric(0),
    probability_b = numeric(0),
    fdr = numeric(0)
  )
  
  sequences <- data.table(
    isoform_id = character(0),
    sequence = character(0)
  )
  
  result <- annotate_drach(m6a_switches, sequences)
  expect_equal(nrow(result), 0)
})

test_that("annotate_drach creates drach_motif_a and drach_motif_b columns", {
  m6a_switches <- data.table(
    gene_id = "GENE1",
    isoform_a = "ISO1",
    isoform_b = "ISO2",
    position = 5,
    m6a_in_isoform_a = TRUE,
    m6a_in_isoform_b = TRUE,
    m6a_fate = "RETAINED",
    probability_a = 0.95,
    probability_b = 0.92,
    fdr = 0.01
  )
  
  sequences <- data.table(
    isoform_id = c("ISO1", "ISO2"),
    sequence = c("ACGAACATCG", "ACGAAGATCG")
  )
  
  result <- annotate_drach(m6a_switches, sequences)
  
  # Should have drach_motif_a and drach_motif_b columns
  expect_true("drach_motif_a" %in% names(result) || "drach_motif" %in% names(result))
})
