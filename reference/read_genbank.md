# Read a GenBank file

Parses a GenBank file and returns a list of sequences, features, and
metadata.

## Usage

``` r
read_genbank(file)
```

## Arguments

- file:

  Path to the GenBank file.

## Value

A list of records, where each record contains:

- metadata:

  A list of metadata with supported fields:

  - `name` (Locus name)

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

- features:

  A data frame of features (key, location, qualifiers)

- sequence:

  The sequence as a string

## Examples

``` r
if (FALSE) { # \dontrun{
  records <- read_genbank("example.gb")
} # }
```
