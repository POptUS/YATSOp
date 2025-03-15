# YATSOp
 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)    [![miss_hit](https://github.com/POptUS/YATSOp/actions/workflows/miss_hit.yml/badge.svg)](https://github.com/POptUS/YATSOp/actions/workflows/miss_hit.yml)

## Yet Another Test Set for Optimization

YATSOp consists of a set of nonlinear optimization problems induced by the solution of systems of (non)linear equations. The defining feature of these problems is that they are of variable input dimension (i.e., the number of decision parameters). This facilitates testing across large ranges of such dimensions.

This page (and perhaps the set) is heavily under construction.

---

## Quickstart

You can find Matlab/Octave versions of the code in [m/](m/) and Python versions of the code in [py/](py/).

In each case, ``calldfofuns`` and ``calldfomidfuns`` provide examples calling particular sets of these problems.

A simple comparison of the two different codebases may be obtained by a Matlab/Octave call to [m/regressiontest.m](m/regressiontest.m) followed by a Python call to [py/compare_fvec.py](py/compare_fvec.py).

## Contributing to YATSOp

Contributions are welcome in a variety of forms; please see [CONTRIBUTING](CONTRIBUTING.rst).

## License 

All code included in YATSOp is open source, with the particular form of license contained in the top-level 
subdirectories.  If such a subdirectory does not contain a LICENSE file, then it is automatically licensed 
as described in the otherwise encompassing YATSOp [LICENSE](/LICENSE).  

## Resources

To seek support or report issues, e-mail:

 * ``poptus@mcs.anl.gov``
