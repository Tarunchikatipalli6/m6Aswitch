test_that("complete workflow: annotate with conditions", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 4),
    ranges = IRanges::IRanges(start = c(100L, 200L, 300L, 400L),
                              end   = c(100L, 200L, 300L, 400L)),
    strand = rep("+", 4))
  gr$transcript_id       <- c("ISO_A", "ISO_A", "ISO_B", "ISO_B")
  gr$transcript_position <- c(10L, 20L, 30L, 40L)
  gr$probability         <- c(0.95, 0.88, 0.92, 0.85)
  gr$condition           <- c("WT", "MUT", "WT", "MUT")

  switches <- data.table::data.table(
    gene_id = "GENE1", isoform_a = "ISO_A", isoform_b = "ISO_B",
    fdr = 0.001, condition_1 = "WT", condition_2 = "MUT", dif = 0.3)

  res <- annotate_m6a_switches_genomic(gr, switches)

  expect_true(data.table::is.data.table(res))
  expect_true(nrow(res) > 0)
  expect_true(all(c("isoform_status", "m6a_fate", "m6a_fate_label",
                    "probability_a", "probability_b") %in% names(res)))
  expect_true(all(res$isoform_status %in%
    c("ISOFORM_A_ONLY", "ISOFORM_B_ONLY", "IN_BOTH_ISOFORMS")))
})

test_that("annotation handles multiple genes", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 4),
    ranges = IRanges::IRanges(start = c(100L, 200L, 500L, 600L),
                              end   = c(100L, 200L, 500L, 600L)),
    strand = rep("+", 4))
  gr$transcript_id       <- c("G1_A", "G1_B", "G2_A", "G2_B")
  gr$transcript_position <- c(10L, 20L, 30L, 40L)
  gr$probability         <- c(0.95, 0.88, 0.92, 0.85)

  switches <- data.table::data.table(
    gene_id     = c("GENE1", "GENE2"),
    isoform_a   = c("G1_A", "G2_A"),
    isoform_b   = c("G1_B", "G2_B"),
    fdr         = c(0.001, 0.01),
    condition_1 = "WT", condition_2 = "MUT", dif = c(0.3, 0.4))

  res <- suppressWarnings(annotate_m6a_switches_genomic(gr, switches))
  expect_setequal(unique(res$gene_id), c("GENE1", "GENE2"))
})

test_that("results are sorted by gene_id then fdr", {
  gr <- GenomicRanges::GRanges(
    seqnames = rep("chr1", 2),
    ranges = IRanges::IRanges(start = c(100L, 200L), end = c(100L, 200L)),
    strand = rep("+", 2))
  gr$transcript_id       <- c("ISO1", "ISO2")
  gr$transcript_position <- c(10L, 20L)
  gr$probability         <- c(0.95, 0.92)

  switches <- data.table::data.table(
    gene_id = "GENE1", isoform_a = "ISO1", isoform_b = "ISO2",
    fdr = 0.01, condition_1 = "WT", condition_2 = "MUT", dif = 0.2)

  res <- suppressWarnings(annotate_m6a_switches_genomic(gr, switches))
  expect_equal(res$gene_id[1], res$gene_id[nrow(res)])
  expect_false(is.unsorted(res$fdr))
})
