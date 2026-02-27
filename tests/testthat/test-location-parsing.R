# Location parsing tests aligned to current rgbio parser behavior

rgbio_parse_location <- function(location) {
  getFromNamespace(".rgbio_parse_feature_location", "rgbio")(location)
}

test_that("location parsing handles simple ranges", {
  result <- rgbio_parse_location("1..100")
  expect_true(is.list(result))
  expect_equal(result$start, 1L)
  expect_equal(result$end, 100L)
  expect_equal(result$strand, "+")
})

test_that("location parsing detects strand correctly", {
  pos_result <- rgbio_parse_location("1..100")
  neg_result <- rgbio_parse_location("complement(1..100)")

  expect_equal(pos_result$strand, "+")
  expect_equal(neg_result$strand, "-")
  expect_false(identical(pos_result, neg_result))
})

test_that("location parsing handles complex locations", {
  join_result <- rgbio_parse_location("join(1..100,200..300)")
  comp_join_result <- rgbio_parse_location("complement(join(1..100,200..300))")
  order_result <- rgbio_parse_location("order(1..100,200..300)")

  expect_equal(join_result$start, 1L)
  expect_equal(join_result$end, 300L)
  expect_equal(comp_join_result$strand, "-")
  expect_equal(order_result$start, 1L)
  expect_equal(order_result$end, 300L)
})

test_that("location parsing from real GenBank files", {
  input_path <- rgbio_biopython_fixture("NC_005816.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_equal(rec$metadata$accession, "NC_005816")
  expect_true("1..9609" %in% rec$features$location)
  expect_true("complement(4815..5888)" %in% rec$features$location)

  parsed_source <- rgbio_parse_location("1..9609")
  expect_equal(parsed_source$start, 1L)
  expect_equal(parsed_source$end, 9609L)
  expect_equal(parsed_source$strand, "+")

  parsed_neg <- rgbio_parse_location("complement(4815..5888)")
  expect_equal(parsed_neg$start, 4815L)
  expect_equal(parsed_neg$end, 5888L)
  expect_equal(parsed_neg$strand, "-")

  if (nrow(rec$features) > 0) {
    for (loc in rec$features$location) {
      expect_true(is.character(loc))
      expect_true(nchar(loc) > 0)
      parsed <- rgbio_parse_location(loc)
      expect_true(is.list(parsed))
      expect_true(is.integer(parsed$start) || is.na(parsed$start))
      expect_true(is.integer(parsed$end) || is.na(parsed$end))
    }
  }
})

test_that("location parsing handles edge case fixture locations", {
  bad_wrap_path <- rgbio_edge_fixture("bad_loc_wrap.gb")
  negative_path <- rgbio_edge_fixture("negative_location.gb")
  one_of_path <- rgbio_biopython_fixture("one_of.gb")

  if (file.exists(bad_wrap_path)) {
    bad <- rgbio_read_wrapper(bad_wrap_path)[[1]]
    expect_equal(nrow(bad$features), 1L)
    parsed <- rgbio_parse_location(bad$features$location[[1]])
    expect_true(is.list(parsed))
    expect_true(is.integer(parsed$start) || is.na(parsed$start))
    expect_true(is.integer(parsed$end) || is.na(parsed$end))
  }

  if (file.exists(negative_path)) {
    neg <- rgbio_read_wrapper(negative_path)[[1]]
    expect_true("complement(-2..492)" %in% neg$features$location)
    parsed <- rgbio_parse_location("complement(-2..492)")
    expect_equal(parsed$start, 2L)
    expect_equal(parsed$end, 492L)
    expect_equal(parsed$strand, "-")
  }

  if (file.exists(one_of_path)) {
    one_of <- rgbio_read_wrapper(one_of_path)[[1]]
    expect_true(any(grepl("one-of\\(1888,1901\\)", one_of$features$location)))
    parsed <- rgbio_parse_location("one-of(1888,1901)")
    expect_equal(parsed$start, 1888L)
    expect_equal(parsed$end, 1901L)
    expect_equal(parsed$strand, "+")
  }
})

test_that("location parsing handles large real feature collections", {
  input_path <- rgbio_biopython_fixture("arab1.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  rec <- records[[1]]

  expect_equal(rec$metadata$accession, "AC007323")
  expect_equal(nchar(rec$sequence), 86436L)
  expect_equal(nrow(rec$features), 19L)

  parsed <- lapply(rec$features$location, rgbio_parse_location)
  expect_true(all(vapply(parsed, function(x) is.integer(x$start) || is.na(x$start), logical(1))))
  expect_true(all(vapply(parsed, function(x) is.integer(x$end) || is.na(x$end), logical(1))))
  expect_true(all(vapply(parsed, function(x) x$strand %in% c("+", "-", "*"), logical(1))))
})

test_that("location parsing validates malformed input safely", {
  invalid_locations <- list("", "invalid", NA_character_, "complement(", "join(,)")

  for (invalid_loc in invalid_locations) {
    parsed <- rgbio_parse_location(invalid_loc)
    expect_true(is.list(parsed))
    expect_true(is.na(parsed$start) || is.integer(parsed$start))
    expect_true(is.na(parsed$end) || is.integer(parsed$end))
  }
})
