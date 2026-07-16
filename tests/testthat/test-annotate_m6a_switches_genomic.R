library(testthat)
library(data.table)
library(GenomicRanges)
library(IRanges)

test_that("annotate_m6a_switches_genomic keeps adjacent loci distinct", {
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(100L, 101L), end = c(100L, 101L)),
    strand = c("+", "+")
  )
  gr$transcript_id <- c("ISO_A", "ISO_B")
  gr$transcript_position <- c(10L, 11L)
  gr$probability <- c(0.9, 0.8)

  switches <- data.table(
    gene_id = "GENE1",
    isoform_a = "ISO_A",
    isoform_b = "ISO_B",
    fdr = 0.01,
    condition_1 = "A",
    condition_2 = "B",
    dif = 0.2
  )

  out <- annotate_m6a_switches_genomic(gr, switches)
  expect_equal(nrow(out), 2)
  expect_setequal(out$m6a_fate, c("LOST", "GAINED"))
})
