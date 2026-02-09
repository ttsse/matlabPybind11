% -------------------------------------------------------------------------
% testMemory.m   -- When running python "out of process" variables passed are 
%                   limited to 2GB. This funciton is used to stop the
%                   program and print an error before a crash is
%                   induced (see RBFsolver, cppFuncTest as examples). 
% Inputs         -- N       -> Positive int. Number of non zeros in A from system 
%                              Ax = b.
%                   M       -> Positive int. Number of elements in b.
%                   tol     -> Positive double. Tolerance used as safety
%                              factor.
% Syntax         -- A = rand(100,100); b = rand(100,1);
%                   testMemory(nnz(A),size(b,1),0.4);      
%
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function [] = testMemory(N,M,tol)

memSize = N*8*1e-9 + M*8*1e-9;

if memSize >= 2 - tol*2 && pyenv().ExecutionMode == "OutOfProcess"
    error("testMemory: System is too large to solve");
end

end