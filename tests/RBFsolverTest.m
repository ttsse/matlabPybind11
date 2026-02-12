% -------------------------------------------------------------------------
% RBFsolverTest.m -- Function-Based Unit tests for RBFsolver function.
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function tests = RBFsolverTest
tests = functiontests(localfunctions);
end
%
% Tests robustness to spurious input (test one input at a time)
%
function testRBFsolverInput(testCase)
    pars = setPars("cppOn",1);
    verifyError(testCase,@() RBFsolver(pars),"RBFsolver:IncorrectType");
    verifyError(testCase,@() RBFsolver,"RBFsolver:IncorrectType");
    verifyError(testCase,@() RBFsolver(pars,pars.cppOn,pars.mode,pars.N),"RBFsolver:IncorrectType");
    verifyError(testCase,@() RBFsolver(pars,pars.cppOn,pars.mode,pars.N,pars.P),"RBFsolver:IncorrectType");
    pars.row = 1;
    verifyError(testCase,@() RBFsolver(pars),"RBFsolver:IncorrectType");
end
%
% Tests output types and fields
%
function testOutputTypes(testCase)
    pars = setPars("cppOn",0,'display',0,'debug',0);
    results = RBFsolver(pars);
    fldNamesRBFsolver = {'u','ue','ucExact','uExact','l2Error','h'}';
    verifyClass(testCase,results,'struct');
    verifyTrue(testCase,all(ismember(fieldnames(results),fldNamesRBFsolver)));
    % Case with debug should include another field (time)
    pars.debug = 1;
    results = RBFsolver(pars);
    fldNamesRBFsolver = {'u','ue','ucExact','uExact','l2Error','h','solveTime'}';
    verifyTrue(testCase,all(ismember(fieldnames(results),fldNamesRBFsolver)));
    for i = 1:length(fldNamesRBFsolver)
        verifyClass(testCase,results.(cell2mat(fldNamesRBFsolver(i))),'double');
    end
end