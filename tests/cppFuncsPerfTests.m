% -------------------------------------------------------------------------
% cppFuncsPerfTests.m -- General performance tests of C++ to matlab connection
%                        and all implemented C++ functions.
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
clearvars;
close all;

pathSet; 
[cpp, np, sp] = importModules; % Import all required modules
N = 1000; q = 2; M = q*N;
Adense = rand(N,N);
Asp = sprand(N,N,0.05,0.1);
b = rand(N,1);
AdenseLS = rand(M,N);
AspLS = sprand(M,N,0.05,0.1);
bLS = rand(M,1);

[BI, BJ, BV] = find(Asp);
row  = np.array(BI - 1).reshape(int32(-1));
col  = np.array(BJ - 1).reshape(int32(-1));
data = np.array(BV,pyargs('dtype', 'float64')).reshape(int32(-1));
pyB = sp.csc_matrix({data,{row, col}}, {int32(size(Asp,1)), int32(size(Asp,2))});

[BI, BJ, BV] = find(AspLS);
row  = np.array(BI - 1).reshape(int32(-1));
col  = np.array(BJ - 1).reshape(int32(-1));
data = np.array(BV,pyargs('dtype', 'float64')).reshape(int32(-1));
pyBLS = sp.csc_matrix({data,{row, col}}, {int32(size(AspLS,1)), int32(size(AspLS,2))});

%% Dense square system solve Matlab
x = Adense\b;

%% Dense square system solve Cpp
x = double(cpp.denseSolve(np.array(Adense),np.array(b),1));

%% Dense least squares system solve Matlab
x = AdenseLS\bLS;

%% Dense least squares system solve Cpp
x = double(cpp.denseSolve(np.array(AdenseLS),np.array(bLS),1));

%% Sparse square system solve Matlab
x = Asp\b;

%% Sparse square system solve C++
x = double(cpp.sparseSolve(pyB,np.array(b,pyargs('dtype', 'float64')),1));

%% Sparse least squares system solve Matlab
x = AspLS\bLS;

%% Sparse least squres system solve C++
x = double(cpp.sparseSolve(pyBLS,np.array(bLS,pyargs('dtype', 'float64')),1));
