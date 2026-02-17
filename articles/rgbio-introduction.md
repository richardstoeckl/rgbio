# Introduction to rgbio

## Overview

`rgbio` provides a high-performance interface to the
[gb-io](https://github.com/moshe/gb-io) Rust crate for reading and
writing GenBank files. It is designed to be fast and memory-efficient
while providing R-friendly data structures.

**Important note:** This package was written primarily by LLMs (“AI”)
under my direction, but it uses the robust Rust `gb-io` crate and is
tested against the NCBI example GenBank files. It is mostly a project to
play around with agents and “AI-tools”, but should provide real value.

## Installation

You can install `rgbio` from GitHub, provided you have the Rust
toolchain installed.

``` r
# install.packages("remotes")
remotes::install_github("richardstoeckl/rgbio")
```

## Basic Usage

### Loading the Package

``` r
library(rgbio)
```

### Writing a GenBank File

To write a GenBank file, you need three components: 1. **Sequence**: A
DNA or RNA string. 2. **Features**: A `data.frame` describing
annotations. 3. **Metadata**: A list of file-level attributes
(definition, accession, etc.).

Let’s create a minimal example sequence.

``` r
# 1. The sequence
seq_dna <- "ATGCGTACGTTAGC"

# 2. Metadata
metadata <- list(
  definition = "Synthetic Example Sequence",
  accession = "EX0001",
  version = "1",
  molecule_type = "DNA",
  topology = "linear",
  division = "SYN",
  date = "01-JAN-2023"
)

# 3. Features
# Note: 'qualifiers' must be a list column where each element is a named character vector.
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

# Preview features
print(features_df)
#>      key location                              qualifiers
#> 1 source    1..14         Synthetic Organism, genomic DNA
#> 2   gene    1..14                             exampleGene
#> 3    CDS    1..14 exampleGene, hypothetical protein, MRTS
```

Now, write it to a temporary file:

``` r
tmp_file <- tempfile(fileext = ".gb")
write_genbank(tmp_file, seq_dna, features_df, metadata)
#> [1] TRUE
```

### Reading a GenBank File

Reading is straightforward. `read_genbank` parses the file and returns a
list of records.

``` r
records <- read_genbank(tmp_file)

# We wrote only one record
length(records)
#> [1] 1

record <- records[[1]]
```

### Inspecting the Data

The returned record has three components matching what we wrote.

**Metadata:**

``` r
str(record$metadata)
#> List of 12
#>  $ name         : chr "EX0001"
#>  $ definition   : chr "Synthetic Example Sequence"
#>  $ accession    : chr "EX0001"
#>  $ version      : chr "1"
#>  $ keywords     : chr(0) 
#>  $ source       : chr ""
#>  $ organism     : chr NA
#>  $ molecule_type: chr "DNA"
#>  $ topology     : chr "linear"
#>  $ division     : chr "SYN"
#>  $ date         : chr "01-JAN-2023"
#>  $ references   : list()
```

**Sequence:**

``` r
record$sequence
#> [1] "ATGCGTACGTTAGC"
```

**Features:**

The features are returned as a tidy `data.frame`.

``` r
print(record$features)
#>      key location   qualifiers
#> 1 source    1..14 Syntheti....
#> 2   gene    1..14  exampleGene
#> 3    CDS    1..14 exampleG....
```

### Supported Metadata Fields

The following metadata fields are supported by
[`write_genbank()`](https://richardstoeckl.github.io/rgbio/reference/write_genbank.md)
and returned by
[`read_genbank()`](https://richardstoeckl.github.io/rgbio/reference/read_genbank.md):

- `name` (Locus name)
- `definition`
- `accession`
- `version`
- `keywords` (character vector)
- `source`
- `organism`
- `molecule_type` (e.g., “DNA”)
- `division`
- `topology` (“linear” or “circular”)
- `date` (format: `DD-MON-YYYY`)
- `references` (list of references; each reference may include
  `description`, `authors`, `consortium`, `title`, `journal`, `pubmed`,
  `remark`)

## Advanced: Complex Locations

`rgbio` supports complex GenBank locations, including joins,
complements, and fuzziness, parsing them into standard location strings.

For advanced manipulations of location strings or arithmetic, you may
need to parse the `location` column further or rely on external
Bioconductor packages. `rgbio` focuses on faithful IO.

## Performance

`rgbio` leverages Rust’s zero-copy parsing where possible and efficient
string handling to outperform pure R implementations, especially for
large multi-record GenBank files.

## Disclaimer

This library is provided under the MIT License. The gb-io Rust crate
package was written by David Leslie and is licensed under the terms of
the MIT License.

This project is in no way affiliated, sponsored, or otherwise endorsed
by the original gb-io authors.
