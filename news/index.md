# Changelog

## rgbio 0.3.0

- Improved GenBank file read performance by ~25%:
  - Improved feature location handling internals with vectorized
    validation and parsing paths.
  - Reduced overhead in
    [`read_gbk()`](https://richardstoeckl.github.io/rgbio/reference/read_gbk.md)
    feature and sequence extraction by reusing vectorized intermediates.
  - Updated write-path validation to validate location vectors in one
    pass for better performance and consistent errors.
  - Expanded internal test coverage for vectorized location validation,
    including scalar/vector error parity checks.

## rgbio 0.2.0

- Initial full release
- Complete rewrite of high-level R functions
  - Added new primary GenBank I/O API:
    [`read_gbk()`](https://richardstoeckl.github.io/rgbio/reference/read_gbk.md)
    and
    [`write_gbk()`](https://richardstoeckl.github.io/rgbio/reference/write_gbk.md).
  - Added `format = "tidy"` and `format = "bioconductor"` output modes
    for reading.
  - Added selective read controls (`sequences`, `features`, `metadata`)
    and record selectors (`records`).
  - Added write support for named character vectors and `DNAStringSet`
    sequence inputs.
  - Added multi-record write orchestration and append-mode writing.
  - Kept
    [`read_genbank()`](https://richardstoeckl.github.io/rgbio/reference/read_genbank.md)
    and
    [`write_genbank()`](https://richardstoeckl.github.io/rgbio/reference/write_genbank.md)
    as deprecated compatibility wrappers.
- Added more extensive testing with test files from gb-io.py, gb-io, and
  biopython

## rgbio 0.1.0

- intermediate release for r-universe

## rgbio 0.0.0.9000

- Initial development version with read/write support for GenBank files.
- R bindings to the Rust `gb-io` crate via extendr.
- Test coverage includes NCBI example GenBank files.
