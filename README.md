# matlabPythonCpp
The present repository implements the Radial Basis Function Partition of Unity method (RBF-PUM) to solve the 
linear elasticity boundary value problem on a thin 3D plank. 

The majority of the code is written in MATLAB with some C++ implementation. The connection is facilitated through
pybind11 bindings as well as the MATLAB Engine API for Python.

## Installation
MATLAB version R2023a is required, as well as the following toolboxes:
- Statistics and Machine Learning Toolbox
- Symbolic Math Toolbox
- Image Processing Toolbox

The conda package manager (or similar) is needed to install the required libraries.
To create the required environment run the following in the terminal:
```
conda create --name TTSSE --file condaEnv.yaml
```
and activate the environment:
```
conda activate TTSSE
```
## Compile the C++ code
To compile the C++ code and create the necessary shared library file to be read by MATLAB use the provided MakeFile:
```
make
```
A file named "cppFuncs.cpython-310-x86_64-linux-gnu.so" should now be available.
To create more C++ functions and add bindings for those functions edit the following [files](/src/cppFuncs)

## Executing the code
Open the MATLAB environment and run main.m
The following parameters can be adjusted as follows:
- Manufactured solution used (man_sol = 'gauss', 'cosine')
- Young's modulus and Poisson's ratio (pars.E > 0, 0 < pars.nu < 0.5)
- Number of patches used (numPatches = [20, 30, 40, 50, 60, 70, 80, 90, 100]); note that any combination and quantity of the given values can be used
- Number of local points used in each patch (nLoc = [56, 84, 120, 165, 220]); note that any combination and quantity of the given values can be used
- Oversampling parameter (pars.oversamp = 3, 5)

The original MATLAB reporitory is not fully provided here. Hence the limitations in parameter definition.
The results from the cover generation algorithm are alternatively given [here](/src/data)

## Results
The code solves the linear elasticity Boundary Value Problem (BVP) for a manufactured solution either in the form of a Gaussian or a cosine function.
The resulting deformation of the body is plotted in MATLAB and should look as follows:


Additional post processing steps can be achieved by running the plotError.m script which plots errors in the discrete $l_2$ norm. 
The resulting error plots should look as follows:

