test_that("read_genbank errors on missing file", {
  missing_path <- tempfile(fileext = ".gb")
  expect_false(file.exists(missing_path))
  expect_error(read_genbank(missing_path), "File not found")
})

test_that("read_genbank errors on empty file", {
  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))
  file.create(tmp)

  expect_error(read_genbank(tmp), "No records found in GenBank file")
})

test_that("read_genbank errors on malformed content", {
  tmp <- tempfile(fileext = ".gb")
  on.exit(unlink(tmp))
  writeLines("INVALID GENBANK CONTENT", tmp)

  expect_error(read_genbank(tmp), "Error reading sequence")
})
