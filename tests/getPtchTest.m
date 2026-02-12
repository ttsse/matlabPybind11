% -------------------------------------------------------------------------
% getPtchTest.m -- Function-Based Unit tests for getPtch function.
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function tests = getPtchTest
tests = functiontests(localfunctions);
end
%
% Run function once to reduce test time (here it is possible to test
% different parameters)
%
function setupOnce(testCase) 
    pars = setPars("mode",'unfitted',"geom",'cube',"dim",2,"del",0.01);
    % pars = setPars("mode",'unfitted',"geom",'cube',"dim",3,"del",0.01);
    % pars = setPars("mode",'unfitted',"geom",'ball',"dim",2,"del",0.01);
    % pars = setPars("mode",'unfitted',"geom",'ball',"dim",3,"del",0.01);
    n = nchoosek(pars.rbfdeg+pars.dim,pars.dim);
    [ptch,dataX,dataY] = getPtch(zeros(1,pars.dim),1,n,pars);
    testCase.TestData.SharedData = {pars,n,ptch,dataX,dataY};
end
%
% Tests robustness to spurious input (test one input at a time)
%
function testGetPtchInput(testCase)
    [pars,n,ptch,dataX,dataY] = testCase.TestData.SharedData{:};
    C = [0,0,5,3];
    verifyError(testCase,@() getPtch(C,1,n,pars),"getPtch:IncorrectType");
    C = zeros(1,pars.dim);
    R = -1;
    verifyError(testCase,@() getPtch(C,R,n,pars),"getPtch:IncorrectType");
    C = zeros(1,pars.dim);
    R = 1;
    n = 0.1;
    verifyError(testCase,@() getPtch(C,R,n,pars),"getPtch:IncorrectType");
    C = zeros(1,pars.dim);
    R = 1;
    n = 5;
    pars.extra = 1;
    verifyError(testCase,@() getPtch(C,R,n,pars),"getPtch:IncorrectType");
    verifyError(testCase,@() getPtch(C,R,n),"getPtch:IncorrectType");
end
%
% Tests output types
%
function testOutputTypes(testCase)
    [pars,n,ptch,dataX,dataY] = testCase.TestData.SharedData{:};
    verifyClass(testCase,ptch,'struct');
    verifyClass(testCase,dataX,'struct');
    verifyClass(testCase,dataY,'struct');
end
%
% Tests output fields and their types
%
function testOutputFields(testCase)
    [pars,n,ptch,dataX,dataY] = testCase.TestData.SharedData{:};
    fldNamesPtchExp = {'C','R','xc','xe'}';
    fldNamesDataXExp = {'nodes','inner','outer','bnd'}';
    fldNamesDataYExp = {'nodes','inner','outer','bnd','Vol','Area'}';
    fldNamesXcExp = {'nodes','globalId'}';
    fldNamesXeExp = {'globalId','nodes'}';
    fldNamesPtch = fieldnames(ptch);
    fldNamesDataX = fieldnames(dataX);
    fldNamesDataY = fieldnames(dataY);
    fldNamesXc = fieldnames(ptch.xc);
    fldNamesXe = fieldnames(ptch.xe);
    verifyTrue(testCase,all(ismember(fldNamesPtchExp,fldNamesPtch)));
    verifyTrue(testCase,all(ismember(fldNamesDataXExp,fldNamesDataX)));
    verifyTrue(testCase,all(ismember(fldNamesDataYExp,fldNamesDataY)));
    verifyTrue(testCase,all(ismember(fldNamesXcExp,fldNamesXc)));
    verifyTrue(testCase,all(ismember(fldNamesXeExp,fldNamesXe)));
    verifyClass(testCase,ptch.xc,'struct');
    verifyClass(testCase,ptch.xe,'struct');
    verifyClass(testCase,ptch.R,'double');
    verifyClass(testCase,ptch.C,'double');
    verifyClass(testCase,dataX.nodes,'double');
    verifyClass(testCase,dataX.inner,'double');
    verifyClass(testCase,dataY.nodes,'double');
    verifyClass(testCase,dataY.inner,'double');
    verifyClass(testCase,dataY.bnd,'double');
    verifyClass(testCase,dataY.Vol,'double');
    verifyClass(testCase,dataY.Area,'double');
    verifyClass(testCase,ptch.xc(1).nodes,'double');
    verifyClass(testCase,ptch.xc(1).globalId,'double');
    verifyClass(testCase,ptch.xe(1).nodes,'double');
    verifyClass(testCase,ptch.xe(1).globalId,'double');
    if ~strcmp(pars.mode,'unfitted')
        verifyClass(testCase,dataY.outer,'double');
        verifyClass(testCase,dataX.bnd,'double');
    else
        verifyEqual(testCase,size(ptch.xc(1).nodes),[n,pars.dim])
        verifyEqual(testCase,size(ptch.xc(1).globalId),[n,1])
    end
end
% 
% Test to see if domain is fully covered
%
function testCover(testCase)
    [pars,n,ptch,dataX,dataY] = testCase.TestData.SharedData{:};
    %
    % Generate many points to see if all are covered
    %
    [dataTest] = getPts(pars.geom,10000,0,zeros(1,pars.dim),1,'fitted',0);
    inCover = zeros(size(dataTest.nodes,1),1);
    for i = 1:size(dataTest.nodes,1)
        distToCentre = sqrt(sum((dataTest.nodes(i,:) - ptch.C).^2,2));
        inCover(i) = any(distToCentre <= ptch.R);
    end
    verifyTrue(testCase,all(inCover));
end