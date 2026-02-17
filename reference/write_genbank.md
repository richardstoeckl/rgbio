# Write a GenBank file

Writes sequences and features to a GenBank format file.

## Usage

``` r
write_genbank(file, sequence, features, metadata = list())
```

## Arguments

- file:

  Path to the output file.

- sequence:

  Character string containing the sequence (DNA/RNA).

- features:

  Data frame of features. Must contain columns:

  - `key`: Feature type (e.g., "CDS", "gene")

  - `location`: Location string (e.g., "1..100")

  - `qualifiers`: List column of named character vectors

- metadata:

  Named list of metadata (optional). Supported fields:

  - `name` (Locus Name)

  - `definition`

  - `accession`

  - `version`

  - `keywords` (character vector)

  - `source`

  - `organism`

  - `molecule_type` (e.g., "DNA")

  - `division`

  - `topology` ("linear" or "circular")

  - `date` (format: `DD-MON-YYYY`)

  - `references` (list of references; each reference may include
    `description`, `authors`, `consortium`, `title`, `journal`,
    `pubmed`, `remark`)

## Value

Logical TRUE on success.

## Examples

``` r
if (FALSE) { # \dontrun{
  meta <- list(definition = "Example Sequence", accession = "AB0001")
  feats <- data.frame(
    key = "source", 
    location = "1..100", 
    qualifiers = I(list(c(organism = "Homo sapiens")))
  )
  write_genbank("out.gb", "ATGC...", feats, meta)
} # }
```
