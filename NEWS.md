# rgbio 0.2.0
- Initial full release
- Complete rewrite of high-level R functions
  - Added new primary GenBank I/O API: `read_gbk()` and `write_gbk()`.
  - Added `format = "tidy"` and `format = "bioconductor"` output modes for reading.
  - Added selective read controls (`sequences`, `features`, `metadata`) and record selectors (`records`).
  - Added write support for named character vectors and `DNAStringSet` sequence inputs.
  - Added multi-record write orchestration and append-mode writing.
  - Kept `read_genbank()` and `write_genbank()` as deprecated compatibility wrappers.
- Added more extensive testing with test files from gb-io.py, gb-io, and biopython

# rgbio 0.1.0

- intermediate release for r-universe

# rgbio 0.0.0.9000

- Initial development version with read/write support for GenBank files.
- R bindings to the Rust `gb-io` crate via extendr.
- Test coverage includes NCBI example GenBank files.
