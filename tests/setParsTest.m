% -------------------------------------------------------------------------
% setParsTest.m -- Function-Based Unit tests for setPars function.
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function tests = setParsTest
tests = functiontests(localfunctions);
end
%
% Tests solution accuracy
%
function testOutput(testCase)
    input = {"dim",3,"memTol",0.1,"q",2,"pdeg",1,"rbfdeg",2,"geom",'cube' ...
             "prob",'Poisson',"method",'FD',"mode",'fitted',"bcMode",'weak', ...
             "scaling",1,"display",1,"mvCentres",0,"psi",'bmp',"phi",'phs',...
             "ep",3,"del",0.4,"P",15,"N",300,"extCoeff",0.4,"cppOn",0,"debug",0};
    pars.dim = 3; pars.geom = 'cube'; pars.prob = 'Poisson'; pars.method = 'FD'; pars.mode = 'fitted';     
    pars.bcMode = 'weak'; pars.scaling = 1; pars.display = 1; pars.mvCentres = 0; pars.psi = 'bmp';                    
    pars.phi = 'phs'; pars.pdeg = 1; pars.ep = 3; pars.memTol = 0.1; pars.del = 0.4; pars.rbfdeg = 2;     
    pars.P = 15; pars.N = 300; pars.extCoeff = 0.4; pars.q = 2; pars.ep = 3; pars.cppOn = 0; pars.debug = 0;
    parsTest = setPars(input{:});
    verifyEqual(testCase,parsTest,pars);
end
%
% Tests robustness to spurious input
%
function testInputsType2(testCase)
    input = {"dim",4};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"memTol",-0.1};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"q",0.9};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"pdeg",-2};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"pdeg",0.5};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"rbfdeg",0};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"geom",'cone'};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"prob",'Heat'};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"method",'FEM'};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"mode",'uCollocation'};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"bcMode",2};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"scaling",1.5};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"display",3};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"mvCentres",'FEM'};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"psi",0.1};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"phi",'w2'};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"ep",-1};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"del",1.1};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"P",12.1};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"N",0};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"extCoeff",0};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"cppOn",0.1};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
    input = {"debug",'yes'};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectType");
end
function testInputsArg(testCase)
    input = {"wrong",4,"phi",'mq',"P",10};
    verifyError(testCase,@() callWithOutput(input),"setPars:IncorrectArg");
end
%
% Make sure if no outputs are provided the function doesnt return warnings
% 
function testNoInput(testCase) 
    input = {"wrong",4,"phi",'mq',"P",10};
    verifyWarningFree(testCase,@() callNoOutput(input));
end
%
% Handle variable outputs
%
function callWithOutput(input)
    [~] = setPars(input{:});
end
function callNoOutput(input)
    setPars(input{:});
end

