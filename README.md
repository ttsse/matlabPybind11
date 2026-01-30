# matlabPythonCpp
The present repository includes the mixed MATLAB and C++ implementation of two RBF based methods (RBF-PUM and RBF-FD). The methods are used to solve a constructed Poisson problem and differential operators are constructed using available MATLAB software. The connection between MATLAB and C++ is facilitated through the pybind11 library and MATLAB Engine API for Python.

The present implementation only works with linux-64 and macOS-64 operating systems given limitations of specific required libraries.

## Installation
MATLAB version R2023a is required, as well as the following toolboxes:
- Statistics and Machine Learning Toolbox
- Symbolic Math Toolbox
- Image Processing Toolbox

The conda package manager is needed to install the required libraries. (See https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). To create the required environment run the following in the terminal:
```
conda create --name TTSSE --file condaEnv.yml
```
and activate the environment:
```
conda activate TTSSE
```
## Compile the C++ code
To compile the C++ code and create the necessary shared library file:
```
make
```
A file named "cppFuncs.cpython-39-x86_64-linux-gnu.so" should now be available.
To create more C++ functions and add bindings for those functions edit the following [files](/src/cppFuncs)

## Download RBF repository
Get repository used to compute derivatives:
```
git clone https://github.com/elisabethl/rbfdiff.git
```

## Tests
The code is tested for both the robustness of the connection between MATLAB and C++ and the accuracy of the solution. The tests are vailable in [files](/tests/) and are as follows:

1. *cppFuncsTest.m*: Verify that results from implemented C++ functions is equal to results from MATLAB functions.
2. *memoryTest.m*: Test to see if inputs to C++ functions can be passed by reference.  
3. *mainInputTest.m*: Ensure that *main.m* is robust to spurious input.
4. *mainErrorTest.m*: Compute numerical error for both methods.
5. *speedTest.m*: Speed comparison between full MATLAB and mixed MATLAB/C++ implementation. 

## Executing the code
Open the MATLAB environment or launch from terminal and run setPaths; followed by main. 
The following parameters can be adjusted in main:

## Results

