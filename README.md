# OOFEM.org
[![Linux](https://github.com/oofem/oofem/actions/workflows/linux.yml/badge.svg)](https://github.com/oofem/oofem/actions/workflows/linux.yml)
[![Windows MSVC](https://github.com/oofem/oofem/actions/workflows/windows.yml/badge.svg)](https://github.com/oofem/oofem/actions/workflows/windows.yml)
[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL%20v2.1-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4339629.svg)](https://doi.org/10.5281/zenodo.4339629)
![GitHub contributors](https://img.shields.io/github/contributors/oofem/oofem) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/oofem/oofem) 
 ![GitHub Downloads (all assets, specific tag)](https://img.shields.io/github/downloads/oofem/oofem/v3.0/total) 



OOFEM is parallel, object-oriented finite element code for solving mechanical, transport and fluid mechanics problems. Read more about its [features here](https://oofem.org/doku.php?id=en:features) or visit [OOFEM web page](https://oofem.org).

OOFEM is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

Copyright (C) 1993 - 2025   Borek Patzak
    
## Getting Started
### What is here
This repository holds source code of the OOFEM.  
```
  OOFEM_TOP_DIR
  |
  |-- doc - contains the "User's guide", sources to generate "Reference manual", 
  |         documents describing the input file specifications, element and
  |         material libraries, and other useful documents. 
  |
  |-- src - source files of all oofem modules
  |   |
  |   |-- core     - sources of the core part of OOFEM
  |   |
  |   |-- sm       - sources of structural analysis module
  |   |
  |   |-- tm       - sources of transport problem module
  |   |
  |   |-- fm       - sources of fluid mechanics module
  |   |
  |   |-- dss      - sources for Direct Sparse Solver
  |   |
  |   |-- main     - contains the sources of main() and implementation of some 
  |                  global functions for oofem, oofeg
  |
  |-- tools   - sources for several utility programs
  |
  |-- tests   - contains tests & benchmarks
  |
  |-- bindings - holds sources to generate OOFEM bindings to Python programming language
```


## First steps
Follow our [Getting Started Guide](http://www.oofem.org/resources/doc/usermanual/index.html) on how to install the OOFEM and make fist steps.

## Documentation & Support
* You may find OOFEM documentation [here](https://www.oofem.org/doku.php?id=en:manual).
* Use [OOFEM forum](https://www.oofem.org/forum/) to post your questions, get support, and much more.

## Authors
See the list of [contributors](http://www.oofem.org/doku.php?id=en:credits) who participated in this project.

## Acknowledgments
http://www.oofem.org/doku.php?id=en:funding




