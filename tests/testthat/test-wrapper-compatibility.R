# Compatibility smoke tests for deprecated wrappers

test_that("read_genbank wrapper remains path-only", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)

  expect_true(is.list(records))
  expect_true(length(records) >= 1)
  expect_true(is.character(records[[1]]$sequence))
  expect_true(is.data.frame(records[[1]]$features))
  expect_true(is.list(records[[1]]$metadata))
})

test_that("read_genbank wrapper rejects non-path input with exact message", {
  expect_error(
    rgbio_read_wrapper(1),
    fixed = TRUE,
    rgbio_expected_messages$file_path_required
  )
})

test_that("write_genbank wrapper writes legacy record", {
  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  features <- data.frame(
    key = "source",
    location = "1..4",
    qualifiers = I(list(c(organism = "synthetic construct"))),
    stringsAsFactors = FALSE
  )

  expect_true(
    rgbio_write_wrapper(
      output_path,
      sequence = "ATGC",
      features = features,
      metadata = list(definition = "Legacy write", accession = "LEGACY001")
    )
  )

  out <- rgbio_read_wrapper(output_path)
  expect_equal(length(out), 1)
  expect_equal(out[[1]]$sequence, "ATGC")
  expect_equal(out[[1]]$metadata$accession, "LEGACY001")
})

test_that("write_genbank wrapper enforces exact legacy validation messages", {
  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  expect_error(
    rgbio_write_wrapper(output_path, sequence = "", features = data.frame(), metadata = list()),
    fixed = TRUE,
    rgbio_expected_messages$sequence_scalar_required
  )

  expect_error(
    rgbio_write_wrapper(output_path, sequence = "ATGC", features = data.frame(), metadata = list()),
    fixed = TRUE,
    rgbio_expected_messages$legacy_features_required
  )
})
