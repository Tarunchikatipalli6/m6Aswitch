test_that("parse_m6anet handles basic input", {
  # Create temporary test data
  test_data <- data.table::data.table(
    transcript_id = c("ENST001", "ENST001", "ENST002"),
    position = c(100, 200, 150),
    probability = c(0.9, 0.7, 0.85)
  )
  
  temp_file <- tempfile(fileext = ".csv")
  data.table::fwrite(test_data, temp_file)
  
  # Parse
  result <- parse_m6anet(temp_file, probability_threshold = 0.5)
  
  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 3)
  expect_true(all(c("transcript_id", "position", "probability") %in% names(result)))
  expect_true(all(result$probability >= 0.5))
  expect_type(result$position, "integer")
  
  unlink(temp_file)
})

test_that("parse_m6anet reports missing required columns before filtering", {
  test_data <- data.table::data.table(
    transcript = "ENST001",
    pos = 100,
    prob = 0.95
  )

  temp_file <- tempfile(fileext = ".csv")
  data.table::fwrite(test_data, temp_file)

  expect_error(
    parse_m6anet(temp_file),
    "must contain columns"
  )

  unlink(temp_file)
})

test_that("parse_m6anet filters by probability threshold", {
  test_data <- data.table::data.table(
    transcript_id = c("ENST001", "ENST001", "ENST002"),
    position = c(100, 200, 150),
    probability = c(0.9, 0.3, 0.85)
  )
  
  temp_file <- tempfile(fileext = ".csv")
  data.table::fwrite(test_data, temp_file)
  
  result <- parse_m6anet(temp_file, probability_threshold = 0.7)
  
  expect_equal(nrow(result), 2)
  expect_true(all(result$probability >= 0.7))
  
  unlink(temp_file)
})
