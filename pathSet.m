%
% Set paths before runing "main.m"
%
addpath("tests/")
addpath(genpath("src/"))
if isfolder("rbfdiff/")
    addpath("rbfdiff")
    addpath("rbfdiff/rbfdir")
    addpath("rbfdiff/rbfqr1D")
    addpath("rbfdiff/rbfqr2D")
    addpath("rbfdiff/rbfqr3D")
    addpath("rbfdiff/utils")
    addpath("rbfdiff/tests")
else
    error("rbfdiff not found, make sure to run: git clone https://github.com/elisabethl/rbfdiff.git in terminal.")
end