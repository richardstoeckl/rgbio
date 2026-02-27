# Read a GenBank file

Reads one or more GenBank records and returns selected components in
either Bioconductor or tidy format.

## Usage

``` r
read_gbk(
  file,
  format = "bioconductor",
  sequences = TRUE,
  features = TRUE,
  metadata = TRUE,
  records = NULL,
  validate = TRUE
)
```

## Arguments

- file:

  Character path to a GenBank file.

- format:

  Output format. One of "bioconductor" or "tidy".

- sequences:

  Logical; include sequence data.

- features:

  Logical; include feature annotations.

- metadata:

  Logical; include record metadata.

- records:

  Integer indices or character accession/name selectors.

- validate:

  Logical; validate parsed records.

## Value

A variable object based on selected components. Returns either a single
object or a named list with `sequences`, `features`, and/or `metadata`.
