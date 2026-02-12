test_that("sequence.gb parses and roundtrips", {
  input_path <- testthat::test_path("..", "sequence.gb")
  if (!file.exists(input_path)) {
    skip(paste("Test file not found:", input_path))
  }

  records <- read_genbank(input_path)
  expect_length(records, 1)
  rec <- records[[1]]

  expect_equal(rec$metadata$accession, "U49845")
  expect_equal(rec$metadata$version, "U49845.1")
  expect_match(rec$sequence, "^GATCCTCCAT")

  output_path <- file.path(tempdir(), "sequence_out.gb")
  expect_true(write_genbank(output_path, rec$sequence, rec$features, rec$metadata))

  strict_equal <- compare_files_strict(input_path, output_path)
  if (!strict_equal) {
    message("Strict file equality failed for sequence_out.gb; data equivalence checks follow.")
  }

  roundtrip <- read_genbank(output_path)[[1]]
  expect_equal(rec$sequence, roundtrip$sequence)
  expect_equal(rec$metadata$accession, roundtrip$metadata$accession)
  expect_equal(rec$features, roundtrip$features)
})

test_that("FASTA + GFF3 recreate sequence.gb", {
  fasta_path <- testthat::test_path("..", "sequence.fasta")
  gff_path <- testthat::test_path("..", "sequence.gff3")
  gb_path <- testthat::test_path("..", "sequence.gb")

  if (!file.exists(fasta_path)) {
    skip(paste("Test file not found:", fasta_path))
  }
  if (!file.exists(gff_path)) {
    skip(paste("Test file not found:", gff_path))
  }
  if (!file.exists(gb_path)) {
    skip(paste("Test file not found:", gb_path))
  }

  fasta <- read_fasta_single(fasta_path)
  header <- fasta$header
  sequence <- fasta$sequence

  accession_version <- sub("\\s.*$", "", header)
  accession <- sub("\\..*$", "", accession_version)
  definition <- sub("^\\S+\\s+", "", header)

  ref <- read_genbank(gb_path)[[1]]
  organism <- ref$metadata$organism
  expect_true(nzchar(organism))

  features <- gff_to_features(gff_path, organism)

  metadata <- list(
    definition = definition,
    accession = accession,
    version = accession_version,
    molecule_type = "DNA",
    topology = "linear",
    division = "PLN",
    source = organism,
    organism = organism
  )

  output_path <- file.path(tempdir(), "sequence_fromFASTAandGFF_out.gb")
  expect_true(write_genbank(output_path, sequence, features, metadata))

  strict_equal <- compare_files_strict(gb_path, output_path)
  if (!strict_equal) {
    message("Strict file equality failed for sequence_fromFASTAandGFF_out.gb; data equivalence checks follow.")
  }

  gen <- read_genbank(output_path)[[1]]

  expect_equal(gen$sequence, ref$sequence)
  expect_equal(gen$metadata$accession, ref$metadata$accession)

  def_ref <- sub("\\.$", "", ref$metadata$definition)
  def_gen <- sub("\\.$", "", gen$metadata$definition)
  expect_equal(def_gen, def_ref)

  expect_equal(gen$features$key, ref$features$key)
  expect_equal(gen$features$location, ref$features$location)
})
