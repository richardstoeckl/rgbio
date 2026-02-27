# Comprehensive loading tests aligned to the read_gbk contract

test_that("read_gbk reads path input in tidy format", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  out <- read_gbk(input_path, format = "tidy")

  expect_true(is.list(out))
  expect_true(all(c("sequences", "features", "metadata") %in% names(out)))
  expect_true(nrow(out$sequences) >= 1)
  expect_true(nrow(out$metadata) >= 1)
  expect_true("sequence" %in% names(out$sequences))
  expect_true("qualifiers" %in% names(out$features))

  expect_equal(out$metadata$accession[[1]], "AY048670")
  expect_equal(out$metadata$version[[1]], "AY048670.1")
  expect_match(out$sequences$sequence[[1]], "^GTCGACTCTA")
  expect_true("source" %in% out$features$type)
  expect_true("CDS" %in% out$features$type)
})

test_that("read_gbk rejects connection-like inputs with exact message", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  con_text <- file(input_path, "r")
  on.exit(close(con_text), add = TRUE)

  expect_error(
    read_gbk(con_text),
    fixed = TRUE,
    rgbio_expected_messages$file_path_required
  )
})

test_that("read_gbk supports record selection by index and accession", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  out_full <- read_gbk(input_path, format = "tidy")
  first_accession <- out_full$metadata$accession[[1]]

  out_by_index <- read_gbk(input_path, format = "tidy", records = 1)
  out_by_accession <- read_gbk(input_path, format = "tidy", records = first_accession)

  expect_equal(nrow(out_by_index$sequences), 1)
  expect_equal(nrow(out_by_accession$sequences), 1)
  expect_equal(out_by_index$sequences$sequence[[1]], out_by_accession$sequences$sequence[[1]])
})

test_that("read_gbk gates compressed file expectations by capability", {
  input_path <- rgbio_basic_fixture("JAOQKG01.1.gb.gz")
  rgbio_require_fixture(input_path)

  if (!rgbio_has_gz_support(input_path)) {
    skip("Gzip read capability not available in current backend/runtime.")
  }

  out <- read_gbk(input_path, format = "tidy")
  expect_true(nrow(out$sequences) >= 1)
  expect_true(all(nchar(out$sequences$sequence) > 0))
})

test_that("read_gbk errors on missing file with exact message", {
  missing_path <- file.path("really", "not", "a", "file", "in", "there.gb")

  expect_error(
    read_gbk(missing_path),
    fixed = TRUE,
    paste0("File not found: ", missing_path)
  )
})

test_that("read_gbk errors on invalid input types with exact message", {
  expect_error(
    read_gbk(1),
    fixed = TRUE,
    rgbio_expected_messages$file_path_required
  )

  expect_error(
    read_gbk(list()),
    fixed = TRUE,
    rgbio_expected_messages$file_path_required
  )

  expect_error(
    read_gbk(NULL),
    fixed = TRUE,
    rgbio_expected_messages$file_path_required
  )
})

test_that("read_gbk enforces records selector contract", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  expect_error(
    read_gbk(input_path, records = list("bad")),
    fixed = TRUE,
    rgbio_expected_messages$records_selector_required
  )
})

test_that("read_gbk enforces at least one selected component", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  expect_error(
    read_gbk(input_path, sequences = FALSE, features = FALSE, metadata = FALSE),
    fixed = TRUE,
    rgbio_expected_messages$nothing_selected
  )
})

test_that("read_gbk propagates parser errors for malformed content", {
  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp), add = TRUE)
  writeLines("LOCUS", tmp)

  expect_error(read_gbk(tmp, format = "tidy"))
})
