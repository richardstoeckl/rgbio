# Write a GenBank file

Deprecated compatibility wrapper for
[`write_gbk()`](https://richardstoeckl.github.io/rgbio/reference/write_gbk.md).

## Usage

``` r
write_genbank(file, sequence, features, metadata = list())
```

## Arguments

- file:

  Path to output file.

- sequence:

  Sequence string.

- features:

  Feature table with columns key, location, qualifiers.

- metadata:

  Metadata list.

## Value

Logical TRUE on success.
