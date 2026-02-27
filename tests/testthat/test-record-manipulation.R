# Record object manipulation tests aligned to R semantics

test_that("record objects follow R copy-on-modify semantics", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec1 <- records[[1]]
  rec2 <- rec1

  original_def <- rec1$metadata$definition
  rec2$metadata$definition <- "Modified definition"

  expect_equal(rec1$metadata$definition, original_def)
  expect_equal(rec2$metadata$definition, "Modified definition")
})

test_that("record objects support deep copying", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec1 <- records[[1]]

  rec2 <- list(
    sequence = rec1$sequence,
    features = if (nrow(rec1$features) > 0) {
      data.frame(rec1$features, stringsAsFactors = FALSE)
    } else {
      data.frame(key = character(), location = character(), qualifiers = I(list()), stringsAsFactors = FALSE)
    },
    metadata = as.list(rec1$metadata)
  )

  original_accession <- rec1$metadata$accession
  rec2$metadata$accession <- "MODIFIED"

  expect_equal(rec1$metadata$accession, original_accession)
  expect_equal(rec2$metadata$accession, "MODIFIED")

  if (nrow(rec1$features) > 0) {
    original_feature_count <- nrow(rec1$features)
    rec2$features <- rec2$features[0, , drop = FALSE]

    expect_equal(nrow(rec1$features), original_feature_count)
    expect_equal(nrow(rec2$features), 0)
  }
})

test_that("record objects maintain structural integrity", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_true(is.character(rec$sequence))
  expect_true(is.data.frame(rec$features))
  expect_true(is.list(rec$metadata))
  expect_true(nchar(rec$sequence) > 0)

  if (nrow(rec$features) > 0) {
    expect_true(all(c("key", "location") %in% names(rec$features)))
  }

  expect_true("accession" %in% names(rec$metadata))
})

test_that("record identity behaves as expected", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec1 <- records[[1]]

  rec2 <- rec1
  expect_identical(rec1, rec2)

  rec3 <- list(
    sequence = rec1$sequence,
    features = rec1$features,
    metadata = rec1$metadata
  )

  expect_equal(rec1$sequence, rec3$sequence)
  expect_false(identical(rec1, rec3))
})

test_that("record objects support modification workflows", {
  input_path <- rgbio_basic_fixture("AY048670.1.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  working_rec <- list(
    sequence = rec$sequence,
    features = if (nrow(rec$features) > 0) {
      data.frame(rec$features, stringsAsFactors = FALSE)
    } else {
      data.frame(key = character(), location = character(), qualifiers = I(list()), stringsAsFactors = FALSE)
    },
    metadata = as.list(rec$metadata)
  )

  working_rec$sequence <- paste0(working_rec$sequence, "AAAA")
  expect_true(nchar(working_rec$sequence) > nchar(rec$sequence))

  new_feature <- data.frame(
    key = "misc_feature",
    location = sprintf("%d..%d", nchar(rec$sequence) + 1, nchar(working_rec$sequence)),
    qualifiers = I(list(c(note = "Added in R"))),
    stringsAsFactors = FALSE
  )

  working_rec$features <- rbind(working_rec$features, new_feature)
  expect_true(nrow(working_rec$features) > nrow(rec$features))

  working_rec$metadata$definition <- paste(working_rec$metadata$definition, "(modified)")
  expect_true(working_rec$metadata$definition != rec$metadata$definition)

  output_path <- tempfile(fileext = ".gb")
  on.exit(unlink(output_path), add = TRUE)

  expect_true(
    rgbio_write_wrapper(output_path, working_rec$sequence, working_rec$features, working_rec$metadata)
  )

  roundtrip <- rgbio_read_wrapper(output_path)
  expect_length(roundtrip, 1)
  expect_equal(roundtrip[[1]]$sequence, working_rec$sequence)
})
