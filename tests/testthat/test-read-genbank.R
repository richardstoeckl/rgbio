test_that("read_genbank errors on missing file", {
  missing_path <- tempfile(fileext = ".gb")
  expect_false(file.exists(missing_path))
  expect_error(rgbio_read_wrapper(missing_path), "File not found")
})

test_that("read_genbank errors on empty file", {
  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))
  file.create(tmp)

  expect_error(rgbio_read_wrapper(tmp), "No records found in GenBank file")
})

test_that("read_genbank errors on malformed content", {
  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))
  writeLines("INVALID GENBANK CONTENT", tmp)

  expect_error(rgbio_read_wrapper(tmp), "Error reading sequence")
})
