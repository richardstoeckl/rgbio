# Write GenBank records

Writes GenBank records from sequence, feature, and metadata components.

## Usage

``` r
write_gbk(
  file,
  sequences,
  features = NULL,
  metadata = NULL,
  append = FALSE,
  validate = TRUE,
  line_width = 80
)
```

## Arguments

- file:

  Character path to output file.

- sequences:

  DNAStringSet or named character vector.

- features:

  GRanges with `type` and `qualifiers` in `mcols()`, or tidy feature
  table with columns `type`, `start`, `end`, `strand`, `qualifiers`.

- metadata:

  DataFrame, data.frame, or list with record metadata.

- append:

  Logical; append to file.

- validate:

  Logical; validate inputs.

- line_width:

  Integer sequence line width.

## Value

Logical TRUE on success.
