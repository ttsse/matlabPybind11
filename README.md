# matlabPybind11
The present repository includes the mixed MATLAB and C++ implementation of two RBF based methods (RBF-PUM and RBF-FD). The methods are used to solve a Poisson problem and differential operators are constructed using available MATLAB software. The connection between MATLAB and C++ is facilitated through the pybind11 library and MATLAB Engine API for Python.

The present implementation works for macOS-64, win-64 and linux-64 operating systems.

## Installation
MATLAB version R2023a is required, as well as the following toolboxes:
- Statistics and Machine Learning Toolbox
- Symbolic Math Toolbox
- Image Processing Toolbox
- Parallel Computing Toolbox

The conda package manager is used to install the required libraries. (See https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). 
## Tools, libraries an relevant lisences
The code used multiple libraries which can be seen in the *.yml files. Note that versions are chosen for cross compatibilit, including with MATLAB. We additionally list all packages and libraries here along with links to their respective lisence files: 
 - gcc 10.4 - [GNU General Public License](https://gcc.gnu.org/onlinedocs/gcc-10.4.0/libstdc++/manual/manual/appendix_gpl.html)
 - gxx 10.4 - [GNU General Public License](https://gcc.gnu.org/onlinedocs/gcc-10.4.0/libstdc++/manual/manual/appendix_gpl.html)
 - make 4.4.1 - [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.html)
 - cmake 4.2.2 - [BSD 3-Clause license](https://gitlab.kitware.com/cmake/cmake/-/blob/master/LICENSE.rst)
 - python 3.9 - [Python Software Foundation & Zero-Clause BSD licenses](https://docs.python.org/3/license.html#bsd0)
 - numpy 1.26 - [BSD 3-Clause License](https://github.com/numpy/numpy/blob/main/LICENSE.txt)
 - Scipy 1.10 - [BSD 3-Clause License](https://github.com/scipy/scipy/blob/main/LICENSE.txt)
 - pybind11 3.0.1 - [BSD 3-Clause License](https://github.com/pybind/pybind11/blob/master/LICENSE)
 - eigen 3.4.0 - [Mozilla Public License 2.0](https://gitlab.com/libeigen/eigen/-/blob/master/LICENSE?ref_type=heads)

Additionally we use the suitesparse 7.10 meta-package. This is a set of packages used for sparse matrix operations and can be found [here](https://github.com/DrTimothyAldenDavis/SuiteSparse). Each package is provided with a different license, all collected [here](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/dev/LICENSE.txt). From the suitesparse repository we specifically use the SuiteSparsQR (SPQR) package [[1]](#1) which also requires the AMD, CAMD [[2]](#2),[[3]](#3), COLAMD, CCOLAMD [[4]](#4),[[5]](#5), CHOLMOD [[6]](#6),[[7]](#7) and the BLAS and LAPACK libraries which are installed separately through OpenBLAS 0.3.31 - [BSD 3-Clause License](https://github.com/OpenMathLib/OpenBLAS/blob/develop/LICENSE).

## Windows
For Windows, Microsoft Visual Studio 2019 is also required, as well as the the Microsoft Visual C++ (MSVC) compiler toolset. (See https://code.visualstudio.com/docs/cpp/config-msvc).
To create the required environment execute the following in the x64 Native Tools Command Prompt for VS:
```
conda create --name TTSSE --file condaEnvWin.yml
```
Next activate the environment:
```
conda activate TTSSE
```
## Linux and macOS
For Linux or macOS execute the following in a terminal:
```
conda create --name TTSSE --file condaEnvLinux.yml
```
Next activate the environment:
```
conda activate TTSSE
```
## Building the project and compiling the shared library
To build the project we use CMake. In the present directory run the following:
```
cmake -S . -B build
cmake --build build
```
The CMake tool will ensure the correct libraries are linked and the correct building and compiling tool are used (based on the OS). 
A file named "cppFuncs.cpython-39-*.so" should now be available.
To create more C++ functions and add bindings for those functions edit the following [files](/src/cppFuncs).

## Download RBF repository
The implementation here uses low level derivative approximation functions from the rbfdiff repository which is licensed under the [GNU General Public License v3.0](https://github.com/elisabethl/rbfdiff/blob/main/LICENSE) and can be cloned:
```
git clone https://github.com/elisabethl/rbfdiff.git
```

## Executing the code
Launch the MATLAB environment from a terminal and run `main`. This executes an example problem with given parameters set by `setPars()`. The example imports and uses the generated C++ library. This code solves the Poisson problem on a 2D circle and should give 4 figures with the points used for the discretisation, the analytical solution, the numerical solution, and the error. All tests and experiments should be executed in the same top/main directory.

## Tests
MATLAB's unit testing framework is used to construct function tests for 3 of the functions used in this implementation, as well as the functions imported from the C++ generated library (see `cppFuncsTest.m`). They are generally tested for robustness to spurious input, expected ouptut and accuracy.  The tests are available in [files](/tests/). All unit tests can be executed from the main directory by first setting the paths through running `pathSet`in the command window and then `runtests("tests")`. If all tests pass, the following should be printed:
```
Totals:
   13 Passed, 0 Failed, 0 Incomplete.
   65.117 seconds testing time.
```

In addition to unit tests there are also script-based performance tests which can be executed by running `pathSet` in the command window and then `runperf("TESTNAME")`. The tests are listed as follows:
 
 - `cppFuncsPerfTests.m` - Testing preformance of implemented C++ functions (Eigen implementation) compared to MATLAB's `mldivide()`.
 - `conGlobMatPerfTests.m` -   Testing preformance of constructing discretisation matrices for both the unfited RBF-FD and RBF-PU methods with and without parallelisation (`parfor`).

Test `cppFuncsTime.m`, is a function based test which outputs wall times for different solve methods (C++ library or `mldivide()`). The function also measures both the solve time prior to any coping of data from MATLAB to C++ and after. This is used to discuss the approximate overhead costs in such operations.

## Experiments
Experiments on performance and accuracy of the implementation can be found [here](src/experiments). These include performance tests of the solver with and without parallelisation (`RBFsolerParExp` which plots wall time against problem size) and convergence studies (`RBFsolverConvExp` which plots error against problem size).

## Armadillo
This branch includes implementation of pybind bindings for sove functions of the Armadillo library licensed under the [Apache License Versoin 2](https://gitlab.com/conradsnicta/armadillo-code/-/blob/15.2.x/LICENSE.txt). To run a simple example follow the steps to create the shared library but while in the armadillo [directory](src\armadillo), then run script `armaEx.m` in matlab.

## References
<a id="1">[1]</a> T. A. Davis, Algorithm 915: SuiteSparseQR: Multifrontal multithreaded rank-revealing sparse QR factorization, ACM Trans. on Mathematical Software, 38(1), 2011, pp. 8:1--8:22. https://doi.org/10.1145/2049662.2049670

<a id="2">[2]</a> P. Amestoy, T. A. Davis, and I. S. Duff, Algorithm 837: An approximate minimum degree ordering algorithm, ACM Trans. on Mathematical Software, 30(3), 2004, pp. 381--388. https://dl.acm.org/doi/abs/10.1145/1024074.1024081

<a id="3">[3]</a> P. Amestoy, T. A. Davis, and I. S. Duff, An approximate minimum degree ordering algorithm, SIAM J. Matrix Analysis and Applications, 17(4), 1996, pp. 886--905. https://doi.org/10.1137/S0895479894278952

<a id="4">[4]</a> T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng, Algorithm 836: COLAMD, an approximate column minimum degree ordering algorithm, ACM Trans. on Mathematical Software, 30(3), 2004, pp. 377--380. https://doi.org/10.1145/1024074.1024080

<a id="5">[5]</a> T. A. Davis, J. R. Gilbert, S. Larimore, E. Ng, A column approximate minimum degree ordering algorithm, ACM Trans. on Mathematical Software, 30(3), 2004, pp. 353--376. https://doi.org/10.1145/1024074.1024079

<a id="6">[6]</a> Y. Chen, T. A. Davis, W. W. Hager, and S. Rajamanickam, Algorithm 887: CHOLMOD, supernodal sparse Cholesky factorization and update/downdate, ACM Trans. on Mathematical Software, 35(3), 2008, pp. 22:1--22:14. https://dl.acm.org/doi/abs/10.1145/1391989.1391995

<a id="7">[7]</a> T. A. Davis and W. W. Hager, Dynamic supernodes in sparse Cholesky update/downdate and triangular solves, ACM Trans. on Mathematical Software, 35(4), 2009, pp. 27:1--27:23. https://doi.org/10.1145/1462173.1462176
