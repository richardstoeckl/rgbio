# rgbio

<!-- badges: start -->
[![R-CMD-check](https://github.com/richardstoeckl/rgbio/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/richardstoeckl/rgbio/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/richardstoeckl/rgbio/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/richardstoeckl/rgbio/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

`rgbio` provides a high-performance interface to the Rust [gb-io](https://github.com/moshe/gb-io) crate for reading and writing GenBank files. It is designed to be fast and memory-efficient while providing R-friendly data structures.

**Important note:** This package was written primarily by LLMs ("AI") under my direction, but it uses the robust Rust `gb-io` crate and is tested against the NCBI example GenBank files. It is mostly a project to play around with agents and "AI-tools", but should provide real value.


## Installation

You can install `rgbio` from GitHub. The Rust toolchain (cargo and rustc) must be available.

```r
# install.packages("remotes")
remotes::install_github("richardstoeckl/rgbio")
```

## Basic Usage

```r
library(rgbio)

# Minimal example
seq_dna <- "ATGCGTACGTTAGC"
metadata <- list(
  definition = "Synthetic Example Sequence",
  accession = "EX0001",
  version = "1",
  molecule_type = "DNA",
  topology = "linear",
  division = "SYN",
  date = "01-JAN-2023"
)

features_df <- data.frame(
  key = c("source", "gene", "CDS"),
  location = c("1..14", "1..14", "1..14"),
  stringsAsFactors = FALSE
)

features_df$qualifiers <- list(
  c(organism = "Synthetic Organism", mol_type = "genomic DNA"),
  c(gene = "exampleGene"),
  c(gene = "exampleGene", product = "hypothetical protein", translation = "MRTS")
)

tmp_file <- tempfile(fileext = ".gb")
write_genbank(tmp_file, seq_dna, features_df, metadata)

records <- read_genbank(tmp_file)
str(records[[1]])
```

## Documentation

The main vignette is available via:

```r
vignette("rgbio-introduction")
```
