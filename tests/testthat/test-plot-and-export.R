library(testthat)
library(data.table)
library(ggplot2)

test_that("delta probability plot keeps LOST/GAINED with not-detected treated as zero", {
  dt <- data.table(
    gene_id = c("G1", "G1", "G1"),
    m6a_fate = c("LOST", "GAINED", "RETAINED"),
    probability_a = c(0.9, NA_real_, 0.4),
    probability_b = c(NA_real_, 0.8, 0.5)
  )

  p <- plot_m6a_switches(dt, plot_type = "delta_prob")
  built <- ggplot_build(p)
  yvals <- built$data[[1]]$y
  expect_equal(length(yvals), 3)
})

test_that("plot_isoform_details filters to requested gene_id", {
  dt <- data.table(
    gene_id = c("G1", "G2"),
    isoform_a = c("A1", "A2"),
    isoform_b = c("B1", "B2"),
    start = c(100L, 200L),
    m6a_fate = c("LOST", "GAINED"),
    m6a_in_isoform_a = c(TRUE, FALSE),
    m6a_in_isoform_b = c(FALSE, TRUE),
    probability_a = c(0.9, NA_real_),
    probability_b = c(NA_real_, 0.9),
    fdr = c(0.01, 0.01)
  )

  p <- plot_isoform_details("G1", dt)
  expect_true(all(p$data$isoform %in% c("A1", "B1")))
})

test_that("export_annotated_switches writes BED9 with optional track line", {
  dt <- data.table(
    gene_id = "G1",
    seqname = "chr1",
    start = 101L,
    end = 101L,
    strand = "+",
    m6a_fate = "LOST",
    probability_a = 0.9,
    probability_b = NA_real_
  )

  out_prefix <- tempfile()
  files <- export_annotated_switches(dt, output_prefix = out_prefix, format = "bed", add_igv_track = TRUE)
  expect_true(file.exists(files$bed))
  bed_lines <- readLines(files$bed)
  expect_match(bed_lines[1], "^track ")
  fields <- strsplit(bed_lines[2], "\t", fixed = TRUE)[[1]]
  expect_equal(length(fields), 9)
})
