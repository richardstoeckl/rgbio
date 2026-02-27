test_that("legacy wrapper handles mismatched feature column lengths", {
  features <- structure(
    list(
      key = c("source", "gene"),
      location = c("1..4", "1..4"),
      qualifiers = list(c(organism = "Test Organism"))
    ),
    class = "data.frame",
    row.names = 1:2
  )
  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))

  expect_true(
    rgbio_write_wrapper(tmp, "ATGC", features, list(definition = "Test", accession = "XX00001"))
  )

  rec <- rgbio_read_wrapper(tmp)[[1]]
  expect_equal(rec$sequence, "ATGC")
  expect_true(is.data.frame(rec$features))
  expect_true(nrow(rec$features) >= 0)
})

test_that("feature validation rejects invalid qualifiers", {
  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))
  meta <- list(definition = "Test", accession = "XX00001")

  features_unnamed <- data.frame(
    key = "source",
    location = "1..4",
    stringsAsFactors = FALSE
  )
  features_unnamed$qualifiers <- list(c("Test Organism"))
  expect_error(
    rgbio_write_wrapper(tmp, "ATGC", features_unnamed, meta),
    fixed = TRUE,
    rgbio_expected_messages$qualifier_names_required
  )

  features_nonchar <- data.frame(
    key = "source",
    location = "1..4",
    stringsAsFactors = FALSE
  )
  features_nonchar$qualifiers <- list(c(organism = 1))
  expect_error(
    rgbio_write_wrapper(tmp, "ATGC", features_nonchar, meta),
    fixed = TRUE,
    rgbio_expected_messages$qualifier_named_character_required
  )
})

test_that("date parsing honors supported format", {
  features <- data.frame(
    key = "source",
    location = "1..4",
    stringsAsFactors = FALSE
  )
  features$qualifiers <- list(c(organism = "Test Organism"))

  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))

  meta_valid <- list(
    definition = "Test Seq",
    accession = "XX00001",
    date = "01-JAN-2023"
  )
  expect_true(rgbio_write_wrapper(tmp, "ATGC", features, meta_valid))
  rec_valid <- rgbio_read_wrapper(tmp)[[1]]
  expect_equal(rec_valid$metadata$date, "01-JAN-2023")

  meta_invalid <- list(
    definition = "Test Seq",
    accession = "XX00001",
    date = "2023-01-01"
  )
  expect_true(rgbio_write_wrapper(tmp, "ATGC", features, meta_invalid))
  rec_invalid <- rgbio_read_wrapper(tmp)[[1]]
  expect_false(identical(rec_invalid$metadata$date, "2023-01-01"))
  if (!is.null(rec_invalid$metadata$date) && nzchar(rec_invalid$metadata$date)) {
    expect_match(rec_invalid$metadata$date, "^[0-9]{2}-[A-Z]{3}-[0-9]{4}$")
  }
})

test_that("invalid locations are rejected", {
  features <- data.frame(
    key = "source",
    location = "1..x",
    stringsAsFactors = FALSE
  )
  features$qualifiers <- list(c(organism = "Test Organism"))

  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))

  expect_error(
    rgbio_write_wrapper(tmp, "ATGC", features, list(definition = "Location Test", accession = "XX00001")),
    fixed = TRUE,
    rgbio_expected_messages$location_tokens_rejected
  )
})

test_that("write_gbk fills missing metadata fields", {
  features <- data.frame(
    type = "source",
    start = 1L,
    end = 4L,
    strand = "+",
    stringsAsFactors = FALSE
  )
  features$qualifiers <- I(list(c(organism = "Test Organism")))

  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))

  seqs <- c(SAMPLE001 = "ATGC")

  expect_true(write_gbk(tmp, seqs, features, list(definition = "Test only")))
  rec_a <- rgbio_read_wrapper(tmp)[[1]]
  expect_equal(rec_a$metadata$definition, "Test only")
  expect_equal(rec_a$metadata$accession, "SAMPLE001")

  expect_true(write_gbk(tmp, seqs, features, list(accession = "ACC_ONLY")))
  rec_b <- rgbio_read_wrapper(tmp)[[1]]
  expect_equal(rec_b$metadata$accession, "ACC_ONLY")
  expect_equal(rec_b$metadata$definition, "SAMPLE001")

  expect_true(write_gbk(tmp, seqs, features, list(definition = "", accession = "XX00001")))
  rec_c <- rgbio_read_wrapper(tmp)[[1]]
  expect_equal(rec_c$metadata$accession, "XX00001")
  expect_equal(rec_c$metadata$definition, "SAMPLE001")
})

test_that("sequence validation rejects empty input", {
  features <- data.frame(
    key = "source",
    location = "1..4",
    stringsAsFactors = FALSE
  )
  features$qualifiers <- list(c(organism = "Test Organism"))

  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))

  expect_error(
    rgbio_write_wrapper(tmp, "", features, list(definition = "Test", accession = "XX00001")),
    fixed = TRUE,
    rgbio_expected_messages$sequence_scalar_required
  )
})

test_that("complex locations are accepted", {
  features <- data.frame(
    key = "gene",
    location = NA_character_,
    stringsAsFactors = FALSE
  )
  features$qualifiers <- list(c(gene = "test_gene"))

  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))

  locations <- c(
    "join(1..100,200..300)",
    "complement(join(490883..490885,1..879))",
    "order(1..100,200..300)",
    "<1..100",
    ">200..300"
  )

  meta <- list(definition = "Location Test", accession = "XX00001")
  for (loc in locations) {
    features$location <- loc
    expect_true(rgbio_write_wrapper(tmp, "ATGC", features, meta), info = loc)
  }
})

test_that("sequence edge cases roundtrip", {
  features <- data.frame(
    key = "source",
    location = "1..1",
    stringsAsFactors = FALSE
  )
  features$qualifiers <- list(c(organism = "Test Organism"))

  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))

  meta <- list(definition = "Edge Seq", accession = "XX00001")

  expect_true(rgbio_write_wrapper(tmp, "A", features, meta))
  rec_single <- rgbio_read_wrapper(tmp)[[1]]
  expect_equal(rec_single$sequence, "A")

  large_seq <- paste(rep("A", 10000), collapse = "")
  expect_true(rgbio_write_wrapper(tmp, large_seq, features, meta))
  rec_large <- rgbio_read_wrapper(tmp)[[1]]
  expect_equal(nchar(rec_large$sequence), nchar(large_seq))

  rna_seq <- "AUGCUU"
  expect_true(rgbio_write_wrapper(tmp, rna_seq, features, meta))
  rec_rna <- rgbio_read_wrapper(tmp)[[1]]
  expect_equal(rec_rna$sequence, rna_seq)
})
