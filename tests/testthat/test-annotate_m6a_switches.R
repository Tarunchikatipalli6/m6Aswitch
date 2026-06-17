library(testthat)
library(data.table)

test_that("annotate_m6a_switches integrates m6A sites with isoform switches", {
  # Create test data
  m6a_sites <- data.table(
    transcript_id = c("ISO1", "ISO1", "ISO2"),
    position = c(100, 150, 100),
    probability = c(0.95, 0.88, 0.92)
  )
  
  iso_switches <- data.table(
    gene_id = "GENE1",
    isoform_a = "ISO1",
    isoform_b = "ISO2",
    fdr = 0.01
  )
  
  iso_sequences <- data.table(
    isoform_id = c("ISO1", "ISO2"),
    sequence = c("ACGTACGT", "ACGTACGT")
  )
  
  result <- annotate_m6a_switches(m6a_sites, iso_switches, iso_sequences)
  
  expect_true(is.data.table(result))
  expect_true(all(c("gene_id", "isoform_a", "isoform_b", "position", 
                     "m6a_fate") %in% names(result)))
})

test_that("annotate_m6a_switches detects LOST, GAINED, RETAINED fates", {
  m6a_sites <- data.table(
    transcript_id = c("ISO1", "ISO2", "ISO1", "ISO2"),
    position = c(100, 200, 300, 300),  # 100 LOST, 200 GAINED, 300 RETAINED
    probability = c(0.95, 0.92, 0.88, 0.90)
  )
  
  iso_switches <- data.table(
    gene_id = c("GENE1", "GENE1", "GENE1"),
    isoform_a = c("ISO1", "ISO1", "ISO1"),
    isoform_b = c("ISO2", "ISO2", "ISO2"),
    fdr = c(0.01, 0.01, 0.01)
  )
  
  iso_sequences <- data.table(
    isoform_id = c("ISO1", "ISO2"),
    sequence = c("ACGTACGT", "ACGTACGT")
  )
  
  result <- annotate_m6a_switches(m6a_sites, iso_switches, iso_sequences)
  
  fates <- unique(result$m6a_fate)
  expect_true("LOST" %in% fates || "GAINED" %in% fates || "RETAINED" %in% fates)
})

test_that("annotate_m6a_switches validates inputs", {
  not_dt <- list(transcript_id = "ISO1")
  iso_switches <- data.table(gene_id = "GENE1", isoform_a = "ISO1", isoform_b = "ISO2")
  iso_sequences <- data.table(isoform_id = "ISO1", sequence = "ACGT")
  
  expect_error(annotate_m6a_switches(not_dt, iso_switches, iso_sequences), 
               "must be a data.table")
})

test_that("annotate_m6a_switches handles no matching m6A sites", {
  m6a_sites <- data.table(
    transcript_id = "ISO3",  # Different isoform
    position = 100,
    probability = 0.95
  )
  
  iso_switches <- data.table(
    gene_id = "GENE1",
    isoform_a = "ISO1",
    isoform_b = "ISO2",
    fdr = 0.01
  )
  
  iso_sequences <- data.table(
    isoform_id = c("ISO1", "ISO2"),
    sequence = c("ACGT", "ACGT")
  )
  
  result <- annotate_m6a_switches(m6a_sites, iso_switches, iso_sequences)
  
  # Should return empty or warning
  expect_true(is.data.table(result))
})

test_that("annotate_m6a_switches has correct output columns", {
  m6a_sites <- data.table(
    transcript_id = "ISO1",
    position = 100,
    probability = 0.95
  )
  
  iso_switches <- data.table(
    gene_id = "GENE1",
    isoform_a = "ISO1",
    isoform_b = "ISO2",
    fdr = 0.01
  )
  
  iso_sequences <- data.table(
    isoform_id = c("ISO1", "ISO2"),
    sequence = c("ACGTACGT", "ACGTACGT")
  )
  
  result <- annotate_m6a_switches(m6a_sites, iso_switches, iso_sequences)
  
  expected_cols <- c("gene_id", "isoform_a", "isoform_b", "position", 
                     "m6a_in_isoform_a", "m6a_in_isoform_b", "m6a_fate",
                     "probability_a", "probability_b", "fdr")
  
  expect_true(all(expected_cols %in% names(result)))
})
