% -------------------------------------------------------------------------
% testMemory.m   -- When running python "out of process" variables passed are 
%                   limited to 2GB. This funciton is used to stop the
%                   program and print an error before a crash is
%                   induced (see RBFsolver, cppFuncTest as examples). 
% Inputs         -- A       -> Double MxN array.
%                   b       -> Double Mx1 array.
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
function [] = checkMemory(A,b,tol)

allVars = whos;
totBytes = 0;
for i = 1:length(allVars)
    totBytes = allVars(i).bytes + totBytes; 
end
totBytes = totBytes/1024^3;

if totBytes >= 2 - tol*2 && pyenv().ExecutionMode == "OutOfProcess"
    error("checkMemory:outOfMemory","checkMemory: System is too large to solve");
end

end