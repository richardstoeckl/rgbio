test_that("all fixture files are referenced by test sources", {
  data_root <- testthat::test_path("..", "data")
  fixture_files <- list.files(data_root, recursive = TRUE, full.names = TRUE)
  fixture_files <- fixture_files[file.info(fixture_files)$isdir %in% FALSE]

  test_sources <- list.files(testthat::test_path(), pattern = "\\.[Rr]$", full.names = TRUE)
  source_lines <- unlist(lapply(test_sources, readLines, warn = FALSE), use.names = FALSE)

  unreferenced <- character(0)
  for (fixture_path in fixture_files) {
    fixture_name <- basename(fixture_path)
    is_referenced <- any(grepl(fixture_name, source_lines, fixed = TRUE))
    if (!is_referenced) {
      unreferenced <- c(unreferenced, fixture_name)
    }
  }

  expect_equal(
    length(unreferenced),
    0L,
    info = paste("Unreferenced fixtures:", paste(sort(unique(unreferenced)), collapse = ", "))
  )
})
