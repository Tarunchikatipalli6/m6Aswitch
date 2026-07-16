library(testthat)
library(data.table)

test_that("annotate_drach validates input", {
  m6a_a <- data.table(transcript_id = "ISO1", position = 10L, probability = 0.9, kmer = "AGACT")
  m6a_b <- copy(m6a_a)
  expect_error(annotate_drach(list(), m6a_a, m6a_b), "must be a data.table")
})

test_that("annotate_drach handles lowercase kmers and does not mutate input by reference", {
  m6a_switches <- data.table(
    gene_id = c("GENE1", "GENE1"),
    isoform_a = c("ISO1", "ISO1"),
    isoform_b = c("ISO2", "ISO2"),
    transcript_position = c(10L, 11L),
    m6a_fate = c("LOST", "GAINED"),
    m6a_in_isoform_a = c(TRUE, FALSE),
    m6a_in_isoform_b = c(FALSE, TRUE),
    probability_a = c(0.95, NA_real_),
    probability_b = c(NA_real_, 0.92),
    fdr = c(0.01, 0.01)
  )
  original <- copy(m6a_switches)

  m6a_condition_a <- data.table(
    transcript_id = "ISO1",
    position = 10L,
    probability = 0.95,
    kmer = "ggact"
  )
  m6a_condition_b <- data.table(
    transcript_id = "ISO2",
    position = 11L,
    probability = 0.92,
    kmer = "AGACT"
  )

  result <- annotate_drach(m6a_switches, m6a_condition_a, m6a_condition_b)

  expect_true("drach_motif" %in% names(result))
  expect_equal(result$drach_motif, c(TRUE, TRUE))
  expect_false("lookup_transcript" %in% names(result))
  expect_false("drach_motif" %in% names(m6a_switches))
  expect_equal(m6a_switches, original)
})
