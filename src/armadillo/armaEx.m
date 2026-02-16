% -------------------------------------------------------------------------
% arnaEx.m -- Armadillo functions (C++ implementation) used to solve linear 
%             systems Ax = b. Including tests of given solution 
%             (wall time / assert conditions).
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
clearvars;
close all;
%
% Parameters: N -> Size of problem to solve, q -> oversampling (for
% rectangular system.
%
N = 1000;
q = 3;
%
% Start pyenv
%
if pyenv().ExecutionMode == "OutOfProcess"
    terminate(pyenv);
end
%
% Import python version from conda environment
%
[~,cmdOut] = system('conda info --base');
cmdOut = strtrim(cmdOut);
if strcmp(computer,"PCWIN64") % GLNXA64 - Linux, MACI64/MACA64 - macOS / silicon
    pyPath = '\envs\TTSSE\python.exe';
else
    pyPath = '/envs/TTSSE/bin/python3';
end
if exist([cmdOut, pyPath])
    pyenv(Version=[cmdOut, pyPath]);
else
    error("cppFuncTest:noConda","cppFuncTest: Install conda (miniconda or anaconda) and generate environment before running this code.")
end
%
% Import python libraries (more can be added here)
%
pyenv("ExecutionMode","OutOfProcess");
armaFuncs = py.importlib.import_module('armaFuncs');
np = py.importlib.import_module('numpy');
%
% Square and rectangular linear problems
%
M = N*q;
A = rand(N,N);
b1 = rand(N,1);
B = rand(M,N);
b2 = rand(M,1);

tic
xSq = double(armaFuncs.solveLSystem(A,b1))';
tSq = toc;

tic
xMat = A\b1;
tMat = toc;

assert(max(abs(xSq-xMat))<1e-10,"error: Solve for square system Ax=b in Armadillo does not match MATLAB solution!");

disp(['armaEx: Armadillo solve time for square system: ', num2str(tSq), 's'])
disp(['armaEx: Matlab solve time for dense system: ', num2str(tMat), 's'])
disp(' ');

tic
xRec = double(armaFuncs.solveLSystem(B,b2))';
tRec = toc;

tic
xMat = B\b2;
tMat = toc;

assert(max(abs(xRec-xMat))<1e-10,"error: Solve for rectangular system Ax=b in Armadillo does not match MATLAB solution!");

disp(['armaEx: Armadillo solve time for rectangular system: ', num2str(tRec), 's'])
disp(['armaEx: Matlab solve time for rectangular system: ', num2str(tMat), 's'])
disp(' ');
%
% Sparse systems (Note: solving sparse under-determined / over-determined systems is currently not supported in Armadillo)
%
Asp = sprand(N,N,0.05,0.1);
tic
xSp = double(armaFuncs.solveSparseSystem(full(Asp),b1))';
tSp = toc;

tic
xMat = Asp\b1;
tMat = toc;

assert(max(abs(xSp-xMat))<1e-10,"error: Solve for sparse square system Ax=b in Armadillo does not match MATLAB solution!");

disp(['armaEx: Armadillo solve time for sparse (square) system: ', num2str(tSp), 's'])
disp(['armaEx: Matlab solve time for sparse (square) system: ', num2str(tMat), 's'])
disp(' ');
