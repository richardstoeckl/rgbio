# Biopython compatibility tests using files from both repositories

rgbio_normalize_legacy_features <- function(features) {
  if (!is.data.frame(features) || !"qualifiers" %in% names(features) || nrow(features) == 0) {
    return(features)
  }

  q_norm <- lapply(features$qualifiers, function(q) {
    if (is.null(q)) {
      return(c(note = ""))
    }
    if (!is.character(q)) {
      q <- as.character(q)
    }
    nm <- names(q)
    if (is.null(nm)) {
      nm <- rep("value", length(q))
    }
    bad <- is.na(nm) | nm == ""
    nm[bad] <- "value"
    names(q) <- nm
    q
  })

  features$qualifiers <- q_norm
  features
}

test_that("handles standard biopython test files", {
  biopython_files <- c("NC_005816.gb", "NT_019265.gb", "arab1.gb", "cor6_6.gb", "iro.gb", "noref.gb", "one_of.gb", "pri1.gb", "protein_refseq.gb")
  existing <- biopython_files[file.exists(vapply(biopython_files, rgbio_biopython_fixture, character(1)))]
  skip_if(length(existing) == 0, "No standard biopython fixtures are available")

  for (filename in existing) {
    records <- rgbio_read_wrapper(rgbio_biopython_fixture(filename))

    if (filename == "cor6_6.gb") {
      expect_equal(length(records), 6L, info = filename)
      expect_equal(
        vapply(records, function(rec) rec$metadata$accession, character(1)),
        c("X55053", "X62281", "M81224", "AJ237582", "L31939", "AF297471"),
        info = filename
      )
      expect_equal(
        vapply(records, function(rec) nchar(rec$sequence), integer(1)),
        c(513L, 880L, 441L, 206L, 282L, 497L),
        info = filename
      )
      next
    }

    rec <- records[[1]]
    expect_equal(length(records), 1L, info = filename)
    expect_true(is.character(rec$sequence), info = filename)
    expect_true(is.data.frame(rec$features), info = filename)
    expect_true(is.list(rec$metadata), info = filename)

    if (filename == "NC_005816.gb") {
      expect_equal(rec$metadata$accession, "NC_005816", info = filename)
      expect_equal(nchar(rec$sequence), 9609L, info = filename)
      expect_equal(nrow(rec$features), 41L, info = filename)
      expect_true("source" %in% rec$features$key, info = filename)
      expect_true("CDS" %in% rec$features$key, info = filename)
      expect_true("complement(4815..5888)" %in% rec$features$location, info = filename)
    }

    if (filename == "NT_019265.gb") {
      expect_equal(rec$metadata$accession, "NT_019265", info = filename)
      expect_equal(nchar(rec$sequence), 0L, info = filename)
      expect_equal(nrow(rec$features), 5L, info = filename)
      expect_true("variation" %in% rec$features$key, info = filename)
    }

    if (filename == "arab1.gb") {
      expect_equal(rec$metadata$accession, "AC007323", info = filename)
      expect_equal(nchar(rec$sequence), 86436L, info = filename)
      expect_equal(nrow(rec$features), 19L, info = filename)
      expect_equal(sort(unique(rec$features$key)), c("CDS", "source"), info = filename)
    }

    if (filename == "iro.gb") {
      expect_equal(rec$metadata$accession, "AL109817", info = filename)
      expect_equal(nchar(rec$sequence), 1326L, info = filename)
      expect_equal(nrow(rec$features), 5L, info = filename)
      expect_true(all(c("source", "gene", "exon", "intron") %in% rec$features$key), info = filename)
    }

    if (filename == "noref.gb") {
      expect_equal(rec$metadata$accession, "NM_006141", info = filename)
      expect_equal(nchar(rec$sequence), 1622L, info = filename)
      expect_equal(nrow(rec$features), 3L, info = filename)
      expect_true(all(c("source", "gene", "CDS") %in% rec$features$key), info = filename)
    }

    if (filename == "one_of.gb") {
      expect_equal(rec$metadata$accession, "U18266", info = filename)
      expect_equal(nchar(rec$sequence), 2509L, info = filename)
      expect_equal(nrow(rec$features), 6L, info = filename)
      expect_true(any(grepl("one-of\\(1888,1901\\)", rec$features$location)), info = filename)
    }

    if (filename == "pri1.gb") {
      expect_equal(rec$metadata$accession, "U05344", info = filename)
      expect_equal(nchar(rec$sequence), 741L, info = filename)
      expect_equal(nrow(rec$features), 5L, info = filename)
      expect_true("repeat_region" %in% rec$features$key, info = filename)
      expect_true("promoter" %in% rec$features$key, info = filename)
    }

    if (filename == "protein_refseq.gb") {
      expect_equal(rec$metadata$accession, "NP_034640", info = filename)
      expect_equal(nchar(rec$sequence), 182L, info = filename)
      expect_equal(nrow(rec$features), 7L, info = filename)
      expect_true("sig_peptide" %in% rec$features$key, info = filename)
      expect_true("mat_peptide" %in% rec$features$key, info = filename)
      expect_true("1..21" %in% rec$features$location, info = filename)
      expect_true("22..182" %in% rec$features$location, info = filename)
    }
  }
})

test_that("handles large genome files", {
  large_files <- c(rgbio_biopython_fixture("NC_000932.gb"), rgbio_basic_fixture("mg1655.gb"))
  existing <- large_files[file.exists(large_files)]
  skip_if(length(existing) == 0, "No large genome fixtures are available")

  for (input_path in existing) {
    records <- rgbio_read_wrapper(input_path)
    rec <- records[[1]]

    expect_length(records, 1)
    expect_true(nchar(rec$sequence) > 10000)
    expect_true(nrow(rec$features) > 10)

    if (grepl("NC_000932\\.gb$", input_path)) {
      expect_equal(rec$metadata$accession, "NC_000932")
      expect_equal(nchar(rec$sequence), 154478L)
      expect_equal(nrow(rec$features), 259L)
    }

    if (grepl("mg1655\\.gb$", input_path)) {
      expect_equal(rec$metadata$accession, "U00096")
      expect_equal(nchar(rec$sequence), 4641652L)
      expect_true(nrow(rec$features) >= 9000)
    }
  }
})

test_that("handles protein GenBank sequences", {
  protein_files <- c(rgbio_biopython_fixture("protein_refseq.gb"), rgbio_biopython_fixture("protein_refseq2.gb"))
  existing <- protein_files[file.exists(protein_files)]
  skip_if(length(existing) == 0, "No protein GenBank fixtures are available")

  for (input_path in existing) {
    records <- rgbio_read_wrapper(input_path)
    rec <- records[[1]]

    expect_length(records, 1)
    expect_equal(rec$metadata$accession, "NP_034640")
    expect_equal(nchar(rec$sequence), 182L)
    expect_true(all(c("source", "Protein", "sig_peptide", "mat_peptide", "CDS") %in% rec$features$key))
    if (grepl("protein_refseq2\\.gb$", input_path)) {
      expect_equal(nrow(rec$features), 9L)
    } else {
      expect_equal(nrow(rec$features), 7L)
    }
  }
})

test_that("handles circular genomes", {
  input_path <- rgbio_basic_fixture("circ.gb")
  rgbio_require_fixture(input_path)

  records <- rgbio_read_wrapper(input_path)
  expect_length(records, 1)

  rec <- records[[1]]
  expect_equal(rec$metadata$accession, "U00096")
  expect_equal(nchar(rec$sequence), 3248L)
  expect_equal(nrow(rec$features), 7L)
  expect_true(grepl("circular", rec$metadata$topology, ignore.case = TRUE))
  expect_true("join(2507..3248,1..2333)" %in% rec$features$location)
})

test_that("handles files with output references", {
  test_pairs <- list(
    c("EU851978.gb", "EU851978_output.gb"),
    c("HM138502.gb", "HM138502_output.gb"),
    c("KF527485.gb", "KF527485_output.gb")
  )

  existing <- Filter(function(pair) {
    file.exists(rgbio_biopython_fixture(pair[1])) && file.exists(rgbio_biopython_fixture(pair[2]))
  }, test_pairs)
  skip_if(length(existing) == 0, "No output-reference fixture pairs are available")

  for (pair in existing) {
    input_records <- rgbio_read_wrapper(rgbio_biopython_fixture(pair[1]))
    output_records <- rgbio_read_wrapper(rgbio_biopython_fixture(pair[2]))
    in_rec <- input_records[[1]]
    out_rec <- output_records[[1]]

    expect_length(input_records, 1)
    expect_length(output_records, 1)
    expect_equal(in_rec$metadata$accession, out_rec$metadata$accession, info = pair[1])
    expect_equal(in_rec$metadata$version, out_rec$metadata$version, info = pair[1])
    expect_equal(in_rec$sequence, out_rec$sequence, info = pair[1])
    expect_equal(nrow(in_rec$features), nrow(out_rec$features), info = pair[1])
    expect_equal(in_rec$features$key, out_rec$features$key, info = pair[1])
    expect_equal(in_rec$features$location, out_rec$features$location, info = pair[1])
  }
})

test_that("handles sequence-only fixtures with stable multi-record output", {
  seq_files <- c(rgbio_basic_fixture("gbvrl1_start.seq"), rgbio_basic_fixture("gbvrl1_start.seq.1"))
  existing <- seq_files[file.exists(seq_files)]
  skip_if(length(existing) == 0, "No sequence-only fixtures are available")

  for (input_path in existing) {
    records <- rgbio_read_wrapper(input_path)
    expect_length(records, 3)
    expect_equal(
      vapply(records, function(rec) rec$metadata$accession, character(1)),
      c("AB000048", "AB000049", "AB000050")
    )
    expect_equal(vapply(records, function(rec) nchar(rec$sequence), integer(1)), c(2007L, 2007L, 1755L))
    expect_equal(vapply(records, function(rec) nrow(rec$features), integer(1)), c(2L, 2L, 2L))
  }
})

test_that("handles GenPept and assembly fixtures with fixture-level expectations", {
  gp_path <- rgbio_biopython_fixture("1MRR_A.gp")
  ds_path <- rgbio_biopython_fixture("DS830848.gb")
  orchid_path <- rgbio_biopython_fixture("ls_orchid.gb")
  pbad_path <- rgbio_biopython_fixture("pBAD30.gb")

  if (file.exists(gp_path)) {
    gp <- rgbio_read_wrapper(gp_path)[[1]]
    expect_equal(gp$metadata$accession, "1MRR_A")
    expect_equal(nchar(gp$sequence), 375L)
    expect_equal(nrow(gp$features), 28L)
    expect_true(all(c("source", "Region", "SecStr", "Site", "Het", "Bond") %in% gp$features$key))
    expect_true("bond(268,272)" %in% gp$features$location)
  }

  if (file.exists(ds_path)) {
    ds <- rgbio_read_wrapper(ds_path)[[1]]
    expect_true(grepl("^DS830848", ds$metadata$accession))
    expect_equal(nchar(ds$sequence), 0L)
    expect_equal(nrow(ds$features), 1L)
    expect_equal(ds$features$key[[1]], "source")
    expect_equal(ds$features$location[[1]], "1..1311")
  }

  if (file.exists(orchid_path)) {
    orchid <- rgbio_read_wrapper(orchid_path)
    expect_length(orchid, 94)
    expect_equal(orchid[[1]]$metadata$accession, "Z78533")
    expect_equal(orchid[[94]]$metadata$accession, "Z78439")
    expect_equal(nchar(orchid[[1]]$sequence), 740L)
    expect_equal(nchar(orchid[[94]]$sequence), 592L)
    expect_true(all(vapply(orchid[c(1, 47, 94)], function(rec) nrow(rec$features), integer(1)) == 5L))
  }

  if (file.exists(pbad_path)) {
    pbad <- rgbio_read_wrapper(pbad_path)[[1]]
    expect_equal(nchar(pbad$sequence), 4923L)
    expect_equal(nrow(pbad$features), 13L)
    expect_true(all(c("rep_origin", "CDS", "promoter", "terminator") %in% pbad$features$key))
    expect_true("complement(1082..1960)" %in% pbad$features$location)
    expect_true("2236..2263" %in% pbad$features$location)
    expect_true("3766..4224" %in% pbad$features$location)
  }
})

test_that("reports parse errors for unsupported TLS/TSA records", {
  tls_path <- rgbio_biopython_fixture("tls_KDHP01000000.gb")
  tsa_path <- rgbio_biopython_fixture("tsa_acropora.gb")

  if (file.exists(tls_path)) {
    expect_error(
      rgbio_read_wrapper(tls_path),
      "TLS|Syntax|MapRes",
      ignore.case = TRUE
    )
  }

  if (file.exists(tsa_path)) {
    expect_error(
      rgbio_read_wrapper(tsa_path),
      "TSA|Syntax|MapRes",
      ignore.case = TRUE
    )
  }
})

test_that("rejects FASTA amino-acid input as non-GenBank", {
  faa_path <- rgbio_basic_fixture("AY048670.1.faa")
  rgbio_require_fixture(faa_path)

  expect_error(
    rgbio_read_wrapper(faa_path),
    "Error reading sequence|No records found|Syntax error",
    ignore.case = TRUE
  )
})

test_that("comprehensive biopython compatibility roundtrip", {
  test_files <- c("arab1.gb", "cor6_6.gb", "iro.gb", "noref.gb", "pri1.gb")
  existing <- test_files[file.exists(vapply(test_files, rgbio_biopython_fixture, character(1)))]
  skip_if(length(existing) == 0, "No roundtrip biopython fixtures are available")

  for (filename in existing) {
    input_path <- rgbio_biopython_fixture(filename)
    original_records <- rgbio_read_wrapper(input_path)
    original_rec <- original_records[[1]]
    features_norm <- rgbio_normalize_legacy_features(original_rec$features)

    output_path <- tempfile(fileext = ".gb")
    on.exit(unlink(output_path), add = TRUE)

    expect_true(
      rgbio_write_wrapper(output_path, original_rec$sequence, features_norm, original_rec$metadata),
      info = filename
    )

    roundtrip_rec <- rgbio_read_wrapper(output_path)[[1]]
    expect_equal(original_rec$sequence, roundtrip_rec$sequence, info = filename)
    if (!is.null(original_rec$metadata$accession)) {
      expect_equal(original_rec$metadata$accession, roundtrip_rec$metadata$accession, info = filename)
    }
    expect_equal(nrow(original_rec$features), nrow(roundtrip_rec$features), info = filename)
  }
})
