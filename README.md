# matlabPythonCpp
The present repository includes the mixed MATLAB and C++ implementation of two RBF based methods (RBF-PUM and RBF-FD). The methods are used to solve a constructed Poisson problem and differential operators are constructed using available MATLAB software. The connection between MATLAB and C++ is facilitated through the pybind11 library and MATLAB Engine API for Python.

The present implementation works for macOS-64, win-64 and linux-64 operating systems.

## Installation
MATLAB version R2023a is required, as well as the following toolboxes:
- Statistics and Machine Learning Toolbox
- Symbolic Math Toolbox
- Image Processing Toolbox

The conda package manager is needed to install the required libraries. (See https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). 
## Libraries used lisences
The following libraries are installed in the conda environment. Note the distribution lisence for each library prior to the use of this code:
 - 
 -

## Windows
For Windows, Microsoft Visual Studio 2019 is required, as well as the the Microsoft Visual C++ (MSVC) compiler toolset. (See https://code.visualstudio.com/docs/cpp/config-msvc).
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
## Building the project and compiling the shared library file
To build the project we have a CMake source file. In this directory run the following:
```
cmake -S . -B build
cmake --build build
```
The CMake tool will ensure the correct libraries are linked and the correct building and compiling tool are used (based on the OS). 
A file named "cppFuncs.cpython-39-*.so" should now be available.
To create more C++ functions and add bindings for those functions edit the following [files](/src/cppFuncs)

## Download RBF repository
The implementation here uses low level derivative approximation functions from the another repository which can be cloned:
```
git clone https://github.com/elisabethl/rbfdiff.git
```

## Tests
The code is tested for both the robustness of the connection between MATLAB and C++ and the accuracy of the solution. The tests are vailable in [files](/tests/) and are as follows:

1. *cppFuncsTest.m*: Verify that results from implemented C++ functions is equal to results from MATLAB functions. Measure the run times.
2. *memoryTest.m*: Test to see if inputs to C++ functions can be passed by reference.  
3. *mainInputTest.m*: Ensure that *main.m* is robust to spurious input.
4. *mainErrorTest.m*: Compute numerical error for both methods.
5. *speedTest.m*: Speed comparison between full MATLAB and mixed MATLAB/C++ implementation. 

## Executing the code
Launch the MATLAB environment from terminal and run setPaths; followed by main. 
The following parameters can be adjusted in main:

## Results

