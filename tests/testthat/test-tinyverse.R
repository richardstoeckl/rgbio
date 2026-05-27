# Tinyverse format and Rust performance parser tests

test_that("read_gbk format = 'base' returns standard R data structures", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  # Read with base format
  res <- read_gbk(input_path, format = "base")

  # Verify structure
  expect_true(is.list(res))
  expect_true(all(c("sequences", "features", "metadata") %in% names(res)))

  # Verify standard data.frame type (no tibbles, no S4 GenomicRanges/DataFrame)
  expect_true(is.data.frame(res$sequences))
  expect_false(inherits(res$sequences, "tbl_df"))
  expect_false(inherits(res$sequences, "DataFrame"))

  expect_true(is.data.frame(res$features))
  expect_false(inherits(res$features, "tbl_df"))
  expect_false(inherits(res$features, "GRanges"))

  expect_true(is.data.frame(res$metadata))
  expect_false(inherits(res$metadata, "tbl_df"))
  expect_false(inherits(res$metadata, "DataFrame"))

  # Verify the column types of features in format = "base"
  expect_type(res$features$start, "integer")
  expect_type(res$features$end, "integer")
  expect_type(res$features$strand, "character")
  expect_true(is.list(res$features$qualifiers))
  expect_true(inherits(res$features$qualifiers, "AsIs"))
})

test_that("Rust location bounds match R regex parser perfectly across complex records", {
  input_path <- rgbio_biopython_fixture("NC_005816.gb")
  rgbio_require_fixture(input_path)

  # Read raw records to compare
  path <- normalizePath(input_path, mustWork = TRUE)
  raw_records <- getFromNamespace(".rgbio_rust_records", "rgbio")(path)

  expect_length(raw_records, 1)
  rec <- raw_records[[1]]

  # Compare Rust parsed bounds/strand against the R-level regex parser
  feat <- rec$features
  r_parsed <- getFromNamespace(".rgbio_parse_feature_location_vec", "rgbio")(feat$location)

  expect_equal(rec$start, r_parsed$start)
  expect_equal(rec$end, r_parsed$end)
  expect_equal(rec$strand, r_parsed$strand)
})

test_that("read_gbk raises helpful errors for missing format packages", {
  # Invalid format validation
  expect_error(
    read_gbk(rgbio_basic_fixture("AY048670.1.gb"), format = "invalid"),
    regexp = "'format' must be one of: 'bioconductor', 'tidy', 'base'"
  )
})
