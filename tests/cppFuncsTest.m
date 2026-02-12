% -------------------------------------------------------------------------
% cppFuncsTest.m -- Function-Based Unit tests to check correctness of functions 
%                   in cppFuncs and testMemory.
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function tests = cppFuncsTest
tests = functiontests(localfunctions);
end
%
% Run import once
%
function setupOnce(testCase) 
    addpath("../")
    [cpp, np, sp] = importModules; % Import all required modules
    N = 1000;
    q = 2;
    M = N*q;
    tol = 1e-9; % Tolerance to compare C++ and Matab solutions (if tests fail try to reduce this)
    testCase.TestData.SharedData = {cpp,np,sp,N,M,tol};
end
%
% Correctness of denseSolve
%
% testMemory(nnz(A),size(b,1),memTol);
function testDenseSolve(testCase)
    [cpp,np,sp,N,M,tol] = testCase.TestData.SharedData{:};
    % Dense square system
    A = rand(N,N);
    b = rand(N,1);
    xDenseCpp = double(cpp.denseSolve(np.array(A),np.array(b)));
    xDenseMat = A\b;
    verifyTrue(testCase,max(abs(xDenseMat-xDenseCpp')) < tol) % Are the solutions "close enough"?
    % Dense least squares system
    B = rand(M,N);
    b = rand(M,1);
    xDenseLScpp = double(cpp.denseSolve(np.array(B),np.array(b)));
    xDenseLSmat = B\b;
    verifyTrue(testCase,max(abs(xDenseLSmat-xDenseLScpp')) < tol) % Are the solutions "close enough"?
end
%
% Correctness of sparseSolve
%
function testSparseSolve(testCase)
    [cpp,np,sp,N,M,tol] = testCase.TestData.SharedData{:};
    % Sparse square system
    A = sprand(N,N,0.05,0.1); % Density of RFB-FD, RBF-PU matrices is approximately 0.05
    b = rand(size(A,1),1);    
    % Making row, col and data arrays 1D. !Necessary for scipy csc_matrix!
    [BI, BJ, BV] = find(A);
    row  = np.array(BI - 1).reshape(int32(-1));
    col  = np.array(BJ - 1).reshape(int32(-1));
    data = np.array(BV,pyargs('dtype', 'float64')).reshape(int32(-1));
    pyB = sp.csc_matrix({data,{row, col}}, {int32(size(A,1)), int32(size(A,2))});
    xSparseCpp = double(cpp.sparseSolve(pyB,np.array(b,pyargs('dtype', 'float64'))));
    xSparseMat = A\b;
    verifyTrue(testCase,max(abs(xSparseMat-xSparseCpp')) < tol) % Are the solutions "close enough"?
    % Sparse square system
    B = sprand(M,N,0.05,0.1); % Density of RFB-FD, RBF-PU matrices is approximately 0.05
    b = rand(size(B,1),1);    
    % Making row, col and data arrays 1D. !Necessary for scipy csc_matrix!
    [BI, BJ, BV] = find(B);
    row  = np.array(BI - 1).reshape(int32(-1));
    col  = np.array(BJ - 1).reshape(int32(-1));
    data = np.array(BV,pyargs('dtype', 'float64')).reshape(int32(-1));
    pyB = sp.csc_matrix({data,{row, col}}, {int32(size(B,1)), int32(size(B,2))});
    xSparseLScpp = double(cpp.sparseSolve(pyB,np.array(b,pyargs('dtype', 'float64'))));
    xSparseLSmat = B\b;
    verifyTrue(testCase,max(abs(xSparseLSmat-xSparseLScpp')) < tol) % Are the solutions "close enough"?
end
%
% Correctness of testMemory (ensuring < 2GB are passed to python -- make sure to stop before python crashes)
%
function testMemoryCheck(testCase)
    [cpp,np] = testCase.TestData.SharedData{:};
    A = rand(13000);
    b = rand(size(A,1),1);
    memTol = 0.4; % Seems to work ok!
    verifyError(testCase,@() checkMemory(A,b,memTol),"checkMemory:outOfMemory")
end
