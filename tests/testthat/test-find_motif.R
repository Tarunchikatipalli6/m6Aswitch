library(testthat)

test_that("find_motif detects DRACH motif", {
  # DRACH: D=A/G/U, R=A/G, A=A, C=C, H=A/C/U
  # Example: GAACA (G=D, A=R, A=A, C=C, A=H)
  sequence <- "GAACA"
  result <- find_motif(sequence, position = 3, context_bp = 2)
  expect_true(result)
})

test_that("find_motif rejects non-DRACH sequence", {
  # CGCCG does not have DRACH pattern
  sequence <- "CGCCG"
  result <- find_motif(sequence, position = 3, context_bp = 2)
  expect_false(result)
})

test_that("find_motif handles edge positions", {
  # Test position at start
  sequence <- "GAACACCG"
  result <- find_motif(sequence, position = 1, context_bp = 2)
  expect_true(is.logical(result))
  
  # Test position at end
  result_end <- find_motif(sequence, position = nchar(sequence), context_bp = 2)
  expect_true(is.logical(result_end))
})

test_that("find_motif validates inputs", {
  expect_error(find_motif("SEQUENCE", position = 100), "within sequence bounds")
  expect_error(find_motif("SEQUENCE", position = 0), "within sequence bounds")
  expect_error(find_motif(c("SEQ1", "SEQ2"), position = 2), "single character string")
})

test_that("find_motif is case insensitive", {
  seq_upper <- "GAACA"
  seq_lower <- "gaaca"
  result_upper <- find_motif(seq_upper, position = 3, context_bp = 2)
  result_lower <- find_motif(seq_lower, position = 3, context_bp = 2)
  expect_equal(result_upper, result_lower)
})

test_that("find_motif with GGACH pattern", {
  # Another valid DRACH: GGACH (G=D, G=R, A=A, C=C, H=any except G)
  sequence <- "GGACA"
  result <- find_motif(sequence, position = 3, context_bp = 2)
  expect_true(result)
})

test_that("find_motif with context_bp parameter", {
  # Test different context windows
  sequence <- "CGAACACCG"
  result_1 <- find_motif(sequence, position = 3, context_bp = 1)
  result_2 <- find_motif(sequence, position = 3, context_bp = 3)
  expect_true(is.logical(result_1))
  expect_true(is.logical(result_2))
})
