# sn_fit_lc

A set of python scripts used to fit supernovae light curves.

```
This software was developed within the LSST DESC using LSST DESC resources, and so meets the criteria 
given in, and is bound by, the LSST DESC Publication Policy for being a "DESC product".
We welcome requests to access code for non-DESC use; if you wish to use the code outside DESC please contact the developers.

```
## Release Status

This code is under development and has not yet been released.

## Feedback, License etc

If you have comments, suggestions or questions, please [write us an issue](https://github.com/LSSTDESC/sn_fit_lc/issues).

This is open source software, available for re-use under the modified BSD license.

```
Copyright (c) 2020, the sn_fit_lc contributors on GitHub, https://github.com/LSSTDESC/sn_fit_lc/graphs/contributors.
All rights reserved.
```

## Content of sn_fit_lc ##
* **docs**: documentation for sphinx
* **\_\_init\_\_.py**
* **input**
* **LICENCE**: licence file
* **README.md**: this readme
* **requirements.txt**: required packages (pip installation) 
* **setup.py**: setup file for pip installation
* [**sn_fit**](doc_package/sn_fit.md): set of scripts to run to fit light curves
* [**sn_fitter**](doc_package/sn_fitter.md): set of fitter to fit the light curves
* **tests**: unit tests
* **version.py**: package version


## Complete tree ##
```bash
|-- LICENCE
|-- README.md
|-- __init__.py
|-- doc_package
|   |-- sn_fit.md
|   |-- sn_fitter.md
|-- docs
|   |-- Makefile
|   |-- api
|   |   |-- sn_fit.mbcov.rst
|   |   |-- sn_fit.process_fit.rst
|   |   |-- sn_fit.rst
|   |   |-- sn_fitter.fit_sn_cosmo.rst
|   |   |-- sn_fitter.fit_sn_fast.rst
|   |   |-- sn_fitter.fit_sncosmo.rst
|   |   |-- sn_fitter.rst
|   |-- conf.py
|   |-- index.rst
|   |-- make.bat
|-- input
|   |-- param.yaml
|-- requirements.txt
|-- setup.py
|-- sn_fit
|   |-- __init__.py
|   |-- mbcov.py
|   |-- process_fit.py
|   |-- version.py
|-- sn_fitter
|   |-- __init__.py
|   |-- fit_sn_cosmo.py
|   `-- fit_sn_fast.py
|-- tests
|   |-- testSNFit.py
|-- version.py

 ```bash
