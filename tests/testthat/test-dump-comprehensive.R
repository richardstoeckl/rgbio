# Comprehensive dumping tests aligned to the write_gbk contract

test_that("write_gbk writes a single record and roundtrips", {
  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  sequence <- c(TEST001 = "ATGCGATCGATCGATCGATCG")
  features <- data.frame(
    type = "CDS",
    start = 1L,
    end = 21L,
    strand = "+",
    qualifiers = I(list(c(gene = "test_gene", product = "test protein"))),
    stringsAsFactors = FALSE
  )
  metadata <- list(
    definition = "Test sequence created in R",
    accession = "TEST001",
    source = "synthetic construct",
    organism = "synthetic construct"
  )

  expect_true(write_gbk(output_path, sequence, features = features, metadata = metadata))

  out <- read_gbk(output_path, format = "tidy")
  expect_equal(nrow(out$sequences), 1)
  expect_equal(out$sequences$sequence[[1]], sequence[[1]])
  expect_equal(out$metadata$accession[[1]], "TEST001")
  expect_true(nrow(out$features) >= 1)
})

test_that("write_gbk appends multiple records with append = TRUE", {
  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  seq1 <- c(REC1 = "ATGCATGC")
  seq2 <- c(REC2 = "GGCCGGCC")

  expect_true(write_gbk(output_path, seq1, features = NULL, metadata = list(accession = "REC1", definition = "Record 1")))
  expect_true(write_gbk(output_path, seq2, features = NULL, metadata = list(accession = "REC2", definition = "Record 2"), append = TRUE))

  out <- read_gbk(output_path, format = "tidy")
  expect_equal(nrow(out$sequences), 2)
  expect_true(all(c("REC1", "REC2") %in% out$metadata$accession))
})

test_that("write_gbk enforces append precondition with exact message", {
  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  expect_error(
    write_gbk(output_path, c(REC1 = "ATGC"), append = TRUE),
    fixed = TRUE,
    rgbio_expected_messages$append_requires_existing
  )
})

test_that("write_gbk errors on invalid file argument with exact message", {
  expect_error(
    write_gbk(NULL, c(REC1 = "ATGC")),
    fixed = TRUE,
    rgbio_expected_messages$file_path_required
  )
})

test_that("write_gbk errors on invalid sequences input with exact message", {
  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  expect_error(
    write_gbk(output_path, list("ATGC")),
    fixed = TRUE,
    rgbio_expected_messages$sequence_set_required
  )
})

test_that("write_gbk errors on invalid features type with exact message", {
  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  expect_error(
    write_gbk(output_path, c(REC1 = "ATGC"), features = "bad"),
    fixed = TRUE,
    rgbio_expected_messages$features_type_required
  )
})

test_that("write_gbk errors on incomplete feature table with exact message", {
  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  bad_features <- data.frame(type = "CDS", stringsAsFactors = FALSE)

  expect_error(
    write_gbk(output_path, c(REC1 = "ATGC"), features = bad_features),
    fixed = TRUE,
    rgbio_expected_messages$features_columns_required
  )
})

test_that("write_gbk errors on invalid metadata type with exact message", {
  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  expect_error(
    write_gbk(output_path, c(REC1 = "ATGC"), metadata = "bad"),
    fixed = TRUE,
    rgbio_expected_messages$metadata_type_required
  )
})

test_that("write_gbk enforces qualifiers structure with exact messages", {
  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  features_unnamed <- data.frame(
    type = "CDS",
    start = 1L,
    end = 4L,
    strand = "+",
    qualifiers = I(list(c("value"))),
    stringsAsFactors = FALSE
  )

  expect_error(
    write_gbk(output_path, c(REC1 = "ATGC"), features = features_unnamed),
    fixed = TRUE,
    rgbio_expected_messages$qualifier_names_required
  )

  features_non_character <- data.frame(
    type = "CDS",
    start = 1L,
    end = 4L,
    strand = "+",
    qualifiers = I(list(c(gene = 1))),
    stringsAsFactors = FALSE
  )

  expect_error(
    write_gbk(output_path, c(REC1 = "ATGC"), features = features_non_character),
    fixed = TRUE,
    rgbio_expected_messages$qualifier_named_character_required
  )
})

test_that("write_gbk warns when non-default line_width is supplied", {
  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  expect_warning(
    write_gbk(output_path, c(REC1 = "ATGC"), line_width = 70),
    "'line_width' is currently ignored by the Rust backend; using backend default.",
    fixed = TRUE
  )
})
