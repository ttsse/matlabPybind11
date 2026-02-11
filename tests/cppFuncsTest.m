% -------------------------------------------------------------------------
% cppFuncsTest.m -- Test all C++ functions and connection with matlab-python.
% Inputs         -- N       -> Positive int. Number of columns used to
%                              generate matrix A.
%                   q       -> Positive int. Oversampling, determines
%                              number of rows in LS problems M = q*N.
%                   memTol  -> Positive double. Tolerance used in testMemory.
%                   display -> Boolean. Display results (including times).
% Outputs        -- tOut    -> Struct. Structure with timings comparing
%                              matlab-python-c++ implementation with pure
%                              matlab implementation.
% Syntax         -- [tOut] = cppFuncsTest(500,2,0.4,1);      
%
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function [tOut] = cppFuncsTest(N,q,memTol,display)
pathSet; 
[cpp, np, sp] = importModules; % Import all required modules
M = N*q;

%
% Solve dense square system
%
A = rand(N,N);
b = rand(N,1);
testMemory(nnz(A),size(b,1),memTol);

tic
xDenseCpp = double(cpp.denseSolve(np.array(A),np.array(b),1));
tCpp = toc;

tic
xDenseMat = A\b;
tMatlab = toc;

if display
    tCppInFunc = xDenseCpp(end);
    xDenseCpp(end) = [];
    disp(['cppFuncTest: Eigen solve time for dense system in function: ', num2str(tCppInFunc), 's'])
    disp(['cppFuncTest: Eigen solve time for dense system: ', num2str(tCpp), 's'])
    disp(['cppFuncTest: Matlab solve time for dense system: ', num2str(tMatlab), 's'])
    disp(' ');
end
tOut.dense = [tCppInFunc, tCpp, tMatlab];
%
% Check if solution from C++ is close enough to Matab solution
%
assert(all(isapprox(xDenseMat,xDenseCpp','loose')),"cppFuncTest: Dense solve solution incorrect!")
%
% Solve dense LS system
%
A = rand(M,N);
b = rand(M,1);
testMemory(nnz(A),size(b,1),memTol);

tic 
xDenseLScpp = double(cpp.denseSolve(np.array(A),np.array(b),1));
tCpp = toc;

tic
xDenseLSmat = A\b;
tMatlab = toc;

if display
    tCppInFunc = xDenseLScpp(end);
    xDenseLScpp(end) = [];
    disp(['cppFuncTest: Eigen solve time for dense LS system in function: ', num2str(tCppInFunc), 's'])
    disp(['cppFuncTest: Eigen solve time for dense LS system: ', num2str(tCpp), 's'])
    disp(['cppFuncTest: Matlab solve time for dense LS system: ', num2str(tMatlab), 's'])
    disp(' ');
end
%
% Check if solution from C++ is close enough to Matab solution
%
assert(all(isapprox(xDenseLSmat,xDenseLScpp','loose')),"cppFuncTest: Dense LS solve solution incorrect!")

tOut.denseLS = [tCppInFunc, tCpp, tMatlab];

%
% Solve sparse system
%
A = sprand(N,N,0.01,0.1);
b = rand(size(A,1),1);
testMemory(nnz(A),size(b,1),memTol);

% Making row, col and data arrays 1D. !Necessary for scipy csc_matrix!
[BI, BJ, BV] = find(A);

row  = np.array(BI - 1).reshape(int32(-1));
col  = np.array(BJ - 1).reshape(int32(-1));
data = np.array(BV,pyargs('dtype', 'float64')).reshape(int32(-1));
pyB = sp.csc_matrix({data,{row, col}}, {int32(size(A,1)), int32(size(A,2))});

tic
xSparseCpp = double(cpp.sparseSolve(pyB,np.array(b,pyargs('dtype', 'float64')),1));
tCpp = toc;

tic
xSparseMat = A\b;
tMatlab = toc;

if display
    tCppInFunc = xSparseCpp(end);
    xSparseCpp(end) = [];
    disp(['cppFuncTest: Eigen solve time for sparse system in function: ', num2str(tCppInFunc), 's'])
    disp(['cppFuncTest: Eigen solve time for sparse system: ', num2str(tCpp), 's'])
    disp(['cppFuncTest: Matlab solve time for sparse system: ', num2str(tMatlab), 's'])
    disp(' ');
end
%
% Check if solution from C++ is close enough to Matab solution
%
assert(all(isapprox(xSparseMat,xSparseCpp','loose')),"cppFuncTest: Sparse solve solution incorrect!")

tOut.sparse = [tCppInFunc, tCpp, tMatlab];

%
% Solve sparse LS system
%
A = sprand(M,N,0.01,0.1);
b = rand(size(A,1),1);
testMemory(nnz(A),size(b,1),memTol);

% Making row, col and data arrays 1D. !Necessary for scipy csc_matrix!
[BI, BJ, BV] = find(A);

row  = np.array(BI - 1).reshape(int32(-1));
col  = np.array(BJ - 1).reshape(int32(-1));
data = np.array(BV,pyargs('dtype', 'float64')).reshape(int32(-1));
pyB = sp.csc_matrix({data,{row, col}}, {int32(size(A,1)), int32(size(A,2))});

tic
xSparseLScpp = double(cpp.sparseSolve(pyB,np.array(b,pyargs('dtype', 'float64')),1));
tCpp = toc;

tic
xSparseLSmat = A\b;
tMatlab = toc;

if display
    tCppInFunc = xSparseLScpp(end);
    xSparseLScpp(end) = [];
    disp(['cppFuncTest: Eigen solve time for sparse LS system in function: ', num2str(tCppInFunc), 's'])
    disp(['cppFuncTest: Eigen solve time for sparse LS system: ', num2str(tCpp), 's'])
    disp(['cppFuncTest: Matlab solve time for sparse LS system: ', num2str(tMatlab), 's'])
    disp(' ');
end
%
% Check if solution from C++ is close enough to Matab solution
%
assert(all(isapprox(xSparseLSmat,xSparseLScpp','loose')),"cppFuncTest: Dense LS solve solution incorrect!")

tOut.sparseLS = [tCppInFunc, tCpp, tMatlab];

end