# The Test Data
## Overview of the tests and their origin

```
tests/
├── testthat/
│   ├── test-load-comprehensive.R          # comprehensive loading tests
│   ├── test-dump-comprehensive.R          # serialization/dumping tests  
│   ├── test-record-manipulation.R         # record object tests
│   ├── test-location-parsing.R            # location parsing tests
│   ├── test-edge-cases.R                  # edge case handling
│   ├── test-biopython-compatibility.R     # Biopython compatibility tests
│   ├── test-read-genbank.R                # rgbio v0.1.0
│   ├── test-internals.R                   # rgbio v0.1.0
│   ├── test-real-data.R                   # rgbio v0.1.0
│   └── helper-io.R                        # rgbio v0.1.0
├── data/
│   ├── basic/
│   │   ├── sequence.gb                    # NCBI Sample GenBank Record (https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
│   │   ├── sequence.fasta                 # NCBI Sample GenBank Record (https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
│   │   ├── sequence.gff3                  # NCBI Sample GenBank Record (https://www.ncbi.nlm.nih.gov/genbank/samplerecord/)
│   │   ├── circ.gb                        # FROM: Rust gb-io/tests/circ.gb
│   │   ├── mg1655.gb                      # FROM: Rust gb-io/tests/mg1655.gb
│   │   ├── AY048670.1.gb                  # FROM: Python gb-io.py/tests/data/AY048670.1.gb
│   │   ├── AY048670.1.faa                 # FROM: Python gb-io.py/tests/data/AY048670.1.faa
│   │   ├── JAOQKG01.1.gb.gz               # FROM: Python gb-io.py/tests/data/JAOQKG01.1.gb.gz
│   │   ├── gbvrl1_start.seq               # FROM: Rust biopython_tests/gbvrl1_start.seq
│   │   └── gbvrl1_start.seq.1             # FROM: Python gb-io.py/tests/data/biopython/gbvrl1_start.seq.1
│   ├── biopython/
│   │   ├── DS830848.gb                    # FROM: Rust biopython_tests/DS830848.gb
│   │   ├── EU851978.gb                    # FROM: Rust biopython_tests/EU851978.gb
│   │   ├── EU851978_output.gb             # FROM: Rust biopython_tests/EU851978_output.gb
│   │   ├── HM138502.gb                    # FROM: Rust biopython_tests/HM138502.gb
│   │   ├── HM138502_output.gb             # FROM: Rust biopython_tests/HM138502_output.gb
│   │   ├── KF527485.gb                    # FROM: Rust biopython_tests/KF527485.gb
│   │   ├── KF527485_output.gb             # FROM: Rust biopython_tests/KF527485_output.gb
│   │   ├── NC_000932.gb                   # FROM: Rust biopython_tests/NC_000932.gb
│   │   ├── NC_005816.gb                   # FROM: Rust biopython_tests/NC_005816.gb
│   │   ├── NT_019265.gb                   # FROM: Rust biopython_tests/NT_019265.gb
│   │   ├── arab1.gb                       # FROM: Rust biopython_tests/arab1.gb
│   │   ├── cor6_6.gb                      # FROM: Rust biopython_tests/cor6_6.gb
│   │   ├── iro.gb                         # FROM: Rust biopython_tests/iro.gb
│   │   ├── ls_orchid.gb                   # FROM: Rust biopython_tests/ls_orchid.gb
│   │   ├── noref.gb                       # FROM: Rust biopython_tests/noref.gb
│   │   ├── one_of.gb                      # FROM: Rust biopython_tests/one_of.gb
│   │   ├── pBAD30.gb                      # FROM: Rust biopython_tests/pBAD30.gb
│   │   ├── pri1.gb                        # FROM: Rust biopython_tests/pri1.gb
│   │   ├── protein_refseq.gb              # FROM: Rust biopython_tests/protein_refseq.gb
│   │   ├── protein_refseq2.gb             # FROM: Rust biopython_tests/protein_refseq2.gb
│   │   ├── 1MRR_A.gp                      # FROM: Python gb-io.py/tests/data/biopython/1MRR_A.gp
│   │   ├── tls_KDHP01000000.gb            # FROM: Python gb-io.py/tests/data/biopython/tls_KDHP01000000.gb
│   │   └── tsa_acropora.gb                # FROM: Python gb-io.py/tests/data/biopython/tsa_acropora.gb
│   ├── edge_cases/
│   │   ├── bad_loc_wrap.gb                # FROM: Rust biopython_tests/bad_loc_wrap.gb
│   │   ├── blank_seq.gb                   # FROM: Rust biopython_tests/blank_seq.gb
│   │   ├── dbsource_wrap.gb               # FROM: Rust biopython_tests/dbsource_wrap.gb
│   │   ├── empty_accession.gb             # FROM: Rust biopython_tests/empty_accession.gb
│   │   ├── empty_feature_qualifier.gb     # FROM: Rust biopython_tests/empty_feature_qualifier.gb
│   │   ├── empty_version.gb               # FROM: Rust biopython_tests/empty_version.gb
│   │   ├── extra_keywords.gb              # FROM: Rust biopython_tests/extra_keywords.gb
│   │   ├── invalid_locus_line_spacing.gb  # FROM: Rust biopython_tests/invalid_locus_line_spacing.gb
│   │   ├── invalid_misc_feature.gb        # FROM: Rust biopython_tests/invalid_misc_feature.gb
│   │   ├── invalid_product.gb             # FROM: Rust biopython_tests/invalid_product.gb
│   │   ├── negative_location.gb           # FROM: Rust biopython_tests/negative_location.gb
│   │   ├── no_end_marker.gb               # FROM: Rust biopython_tests/no_end_marker.gb
│   │   ├── origin_line.gb                 # FROM: Rust biopython_tests/origin_line.gb
│   │   └── wrong_sequence_indent.gb       # FROM: Rust biopython_tests/wrong_sequence_indent.gb
└── testthat.R                             # rgbio v0.1.0
```

## Test Data Licensing and Attribution

### Sources

The test files in this directory are derived from multiple sources:

1. **Biopython test files** (in `biopython/` and `edge_cases/` directories): Originally from the Biopython project test suite, copied via:
   - Rust `gb-io` repository: https://github.com/dlesl/gb-io/tree/master/tests/biopython_tests/
   - Python `gb-io.py` repository: https://github.com/althonos/gb-io.py/tree/main/tests/data/

2. **Additional test files** (in `basic/` directory): From the Rust `gb-io` and Python `gb-io.py` repositories' test suites

3. **NCBI Sample GenBank Record** (in `basic/` directory): Annotated sample GenBank record (accession number U49845) in its GenBank Flat File format (https://www.ncbi.nlm.nih.gov/genbank/samplerecord/).

### Copyright and Licensing

#### Biopython Test Files

The Biopython test files are licensed under the **Biopython License Agreement** unless otherwise specified in individual file headers. Some files are dual-licensed under either the "Biopython License Agreement" or the "BSD 3-Clause License".

**Copyright Holders** (as noted in the Python gb-io.py test suite):
- Brad Chapman (2001-2004)
- Peter Cock (2007-2016, 2015-2016) 
- Kai Blin (2013)
- Sergio Valqui (2019)
- The Biopython Contributors (1999-2025)

#### Biopython License Agreement

Permission to use, copy, modify, and distribute this software and its documentation with or without modifications and for any purpose and without fee is hereby granted, provided that any copyright notices appear in all copies and that both those copyright notices and this permission notice appear in supporting documentation, and that the names of the contributors or copyright holders not be used in advertising or publicity pertaining to distribution of the software without specific prior permission.

THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

#### BSD 3-Clause License (Alternative License)

Copyright (c) 1999-2025, The Biopython Contributors
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

### Usage in This Package

These test files are used solely for testing the `rgbio` R package's GenBank parsing functionality. They are included here under the terms of the above licenses to ensure comprehensive testing compatibility with the broader bioinformatics ecosystem.

### Acknowledgments

We gratefully acknowledge:
- The Biopython project and its contributors for creating and maintaining these valuable test datasets
- The maintainers of the Rust `gb-io` and Python `gb-io.py` packages for organizing and curating these test files
- All the original copyright holders listed above for their contributions to the bioinformatics community

