# Edge case tests based on Rust gb-io biopython_tests files

test_that("handles blank sequence files", {
  input_path <- rgbio_edge_fixture("blank_seq.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  expect_length(records, 1)
  rec <- records[[1]]

  expect_equal(rec$metadata$accession, "NP_001832")
  expect_equal(nchar(rec$sequence), 360L)
  expect_equal(nrow(rec$features), 4L)
  expect_true(all(c("source", "Protein", "Region", "CDS") %in% rec$features$key))
})

test_that("handles files with bad location wrapping", {
  input_path <- rgbio_edge_fixture("bad_loc_wrap.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_length(records, 1)
  expect_equal(nchar(rec$sequence), 6000L)
  expect_equal(nrow(rec$features), 1L)
  expect_true("CDS" %in% rec$features$key)
})

test_that("handles files with dbsource wrapping", {
  input_path <- rgbio_edge_fixture("dbsource_wrap.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  expect_length(records, 1)
  rec <- records[[1]]

  expect_equal(rec$metadata$accession, "P01485")
  expect_equal(nchar(rec$sequence), 64L)
  expect_equal(nrow(rec$features), 7L)
  expect_equal(sum(rec$features$key == "Bond"), 4L)
  expect_equal(sum(rec$features$key == "Site"), 1L)
  expect_true("bond(12,63)" %in% rec$features$location)
})

test_that("handles files with empty accession", {
  input_path <- rgbio_edge_fixture("empty_accession.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]
  accession <- rec$metadata$accession

  expect_equal(nchar(rec$sequence), 6497L)
  expect_equal(nrow(rec$features), 1L)
  expect_true("source" %in% rec$features$key)
  expect_true(is.null(accession) || is.na(accession) || nchar(accession) == 0)
})

test_that("handles files with empty feature qualifiers", {
  input_path <- rgbio_edge_fixture("empty_feature_qualifier.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_equal(rec$metadata$accession, "AB070938")
  expect_equal(nchar(rec$sequence), 6497L)
  expect_equal(nrow(rec$features), 1L)
  expect_true("qualifiers" %in% names(rec$features))
  expect_equal(rec$features$key[[1]], "source")
})

test_that("handles files with empty version", {
  input_path <- rgbio_edge_fixture("empty_version.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]
  version <- rec$metadata$version

  expect_equal(nchar(rec$sequence), 6497L)
  expect_equal(nrow(rec$features), 1L)
  expect_true("source" %in% rec$features$key)
  expect_true(is.null(version) || is.na(version) || nchar(version) == 0)
})

test_that("handles files with extra keywords", {
  input_path <- rgbio_edge_fixture("extra_keywords.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_equal(rec$metadata$accession, "AL138972")
  expect_equal(nchar(rec$sequence), 154329L)
  expect_equal(nrow(rec$features), 25L)
  expect_true(sum(rec$features$key == "gene") >= 1)
  expect_true(sum(rec$features$key == "CDS") >= 1)
})

test_that("handles files with invalid locus line spacing", {
  input_path <- rgbio_edge_fixture("invalid_locus_line_spacing.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_length(records, 1)
  expect_equal(rec$metadata$accession, "AB070938")
  expect_equal(nchar(rec$sequence), 6497L)
  expect_equal(nrow(rec$features), 1L)
})

test_that("handles files with invalid misc features", {
  input_path <- rgbio_edge_fixture("invalid_misc_feature.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_length(records, 1)
  expect_equal(rec$metadata$accession, "AB070938")
  expect_equal(nchar(rec$sequence), 6497L)
  expect_equal(nrow(rec$features), 2L)
  expect_true(any(grepl("^misc_feature", rec$features$key)))
})

test_that("handles files with invalid product", {
  input_path <- rgbio_edge_fixture("invalid_product.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_length(records, 1)
  expect_equal(rec$metadata$accession, "AB070938")
  expect_equal(nchar(rec$sequence), 6497L)
  expect_equal(nrow(rec$features), 3L)
  expect_true("contain" %in% rec$features$key)
  expect_true("CDS" %in% rec$features$key)
})

test_that("handles files with negative locations", {
  input_path <- rgbio_edge_fixture("negative_location.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_length(records, 1)
  expect_equal(rec$metadata$accession, "AB070938")
  expect_equal(nchar(rec$sequence), 6497L)
  expect_equal(nrow(rec$features), 2L)
  expect_true("complement(-2..492)" %in% rec$features$location)
})

test_that("handles files with no end marker", {
  input_path <- rgbio_edge_fixture("no_end_marker.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_length(records, 1)
  expect_equal(rec$metadata$accession, "AB070938")
  expect_equal(nchar(rec$sequence), 6497L)
  expect_equal(nrow(rec$features), 1L)
})

test_that("handles files with origin line issues", {
  input_path <- rgbio_edge_fixture("origin_line.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  expect_length(records, 1)
  rec <- records[[1]]

  expect_equal(rec$metadata$accession, "NC_002678")
  expect_equal(nchar(rec$sequence), 180L)
  expect_equal(nrow(rec$features), 2L)
  expect_true("20..120" %in% rec$features$location)
})

test_that("handles files with wrong sequence indentation", {
  input_path <- rgbio_edge_fixture("wrong_sequence_indent.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_length(records, 1)
  expect_equal(rec$metadata$accession, "AB070938")
  expect_equal(nchar(rec$sequence), 6497L)
  expect_equal(nrow(rec$features), 1L)
})

test_that("roundtrip testing with edge case files", {
  edge_case_files <- c(
    "empty_accession.gb",
    "empty_feature_qualifier.gb",
    "empty_version.gb",
    "origin_line.gb"
  )

  existing <- edge_case_files[file.exists(vapply(edge_case_files, rgbio_edge_fixture, character(1)))]
  skip_if(length(existing) == 0, "No roundtrip edge-case fixtures are available")

  for (filename in existing) {
    input_path <- rgbio_edge_fixture(filename)

    records <- rgbio_read_wrapper(input_path)
    expect_equal(length(records), 1L, info = filename)

    rec <- records[[1]]
    output_path <- tempfile(fileext = ".gb")
    on.exit(unlink(output_path), add = TRUE)

    expect_true(
      rgbio_write_wrapper(output_path, rec$sequence, rec$features, rec$metadata),
      info = filename
    )

    roundtrip <- rgbio_read_wrapper(output_path)
    expect_equal(length(roundtrip), 1L, info = filename)
    expect_equal(roundtrip[[1]]$sequence, rec$sequence, info = filename)
    expect_equal(nrow(roundtrip[[1]]$features), nrow(rec$features), info = filename)

    if (!is.null(rec$metadata$accession) && !is.na(rec$metadata$accession) && nzchar(rec$metadata$accession)) {
      expect_equal(roundtrip[[1]]$metadata$accession, rec$metadata$accession, info = filename)
    }
  }
})
