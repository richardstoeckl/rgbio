# rgbio

<!-- badges: start -->
[![R-CMD-check](https://github.com/richardstoeckl/rgbio/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/richardstoeckl/rgbio/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/richardstoeckl/rgbio/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/richardstoeckl/rgbio/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

`rgbio` provides a high-performance interface to the Rust [gb-io](https://github.com/moshe/gb-io) crate for reading and writing GenBank files. It is designed to be fast and memory-efficient while providing R-friendly data structures.

**Important note:** This package was written primarily by LLMs ("AI") under my direction, but it uses the robust Rust `gb-io` crate and is tested against ~50 diverse GenBank files with many edge cases. 
It is a project for me to play around with agentic coding, but provides real value as it is one of the only ways to write GenBank files in R.


## Installation

The `rgbio` package is not available on CRAN (for now), because it depends on a Rust crate. 
You can install it from the R-universe repository without having installed Rust or any Rust toolchain, as there are binary versions available for Windows, macOS, and Linux.

```
install.packages("rgbio", 
                 repos = c("https://richardstoeckl.r-universe.dev", 
                           "https://cloud.r-project.org"))
```

If there is no pre-built binary available for your system, or you want the latest development version, you can install `rgbio` from GitHub, provided you have the Rust toolchain installed.
You can find information on how to install Rust at https://github.com/r-rust/hellorust.
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
  type = c("source", "gene", "CDS"),
  start = c(1L, 1L, 1L),
  end = c(14L, 14L, 14L),
  strand = c("+", "+", "+"),
  stringsAsFactors = FALSE
)

features_df$qualifiers <- list(
  c(organism = "Synthetic Organism", mol_type = "genomic DNA"),
  c(gene = "exampleGene"),
  c(gene = "exampleGene", product = "hypothetical protein", translation = "MRTS")
)

tmp_file <- tempfile(fileext = ".gb")
write_gbk(
  file = tmp_file,
  sequences = c(EX0001 = seq_dna),
  features = features_df,
  metadata = metadata
)

records <- read_gbk(tmp_file, format = "tidy")
str(records)
```

## What Is Required vs Optional for `write_gbk()`

Minimum required inputs:

- `file`: output file path
- `sequences`: non-empty named character vector or `DNAStringSet`

Optional inputs:

- `features = NULL`
- `metadata = NULL`
- `append = FALSE`
- `validate = TRUE`
- `line_width = 80`

If `metadata` is omitted, `rgbio` fills defaults per record:

- `name`, `definition`, `accession` from the sequence name
- `molecule_type = "DNA"`

`append = TRUE` requires that the target file already exists and is a valid GenBank file.

## Minimal Tidy Example

```r
library(rgbio)

out <- tempfile(fileext = ".gb")

write_gbk(
  file = out,
  sequences = c(min1 = "ATGC")
)

read_gbk(out, format = "tidy")
```

## Minimal Bioconductor Example

```r
library(rgbio)

out <- tempfile(fileext = ".gb")

seqs <- Biostrings::DNAStringSet(c(min2 = "ATGCGG"))
write_gbk(
  file = out,
  sequences = seqs
)

read_gbk(out, format = "bioconductor")
```

## Documentation

The main vignette is available via:

```r
vignette("rgbio-introduction")
```
