% -------------------------------------------------------------------------
% conGlobMatPerfTests.m -- General script based performance tests for parallelisation of
%                          conGlobMat. Run with runperf("conGlobMatPerfTests");
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
clearvars;
close all;
%
% Setup RBF-FD parameters
%
parsFD = setPars("geom",'ball',"N",1000,"display",0,"debug",0,"rbfdeg",3,"method","FD");
n = nchoosek(parsFD.rbfdeg+parsFD.dim,parsFD.dim);
C = zeros(1,parsFD.dim);
R = 1;
[dataXFD] = getPts(parsFD.geom,parsFD.N,n,C,R,parsFD.mode,parsFD.extCoeff);
xc = dataXFD.nodes;         
parsFD.M = ceil(parsFD.q*parsFD.N);   
[dataYFD] = getPts(parsFD.geom,parsFD.M,n,C,R,"fitted",0);
dataYFD = movePts(dataXFD,dataYFD);
xe = dataYFD.nodes;
xc = xc(unique(knnsearch(xc,xe,'K',ceil(parsFD.extCoeff*n+parsFD.ep))),:);
parsFD.N = size(xc,1);
dataXFD.nodes = xc;
ptchFD = knnsearch(xc,xe,'K',1); % This is a list of closest centre points to eval points (centres)
parsFD.n = n;
parsFD.M = size(dataYFD.nodes,1);
parsFD.N = size(dataXFD.nodes,1); 
%
% Setup RBF-PUM parameters
%
parsPUM = setPars("geom",'ball',"P",100,"display",0,"debug",0,"rbfdeg",3,"method",'PUM');
n = nchoosek(parsPUM.rbfdeg+parsPUM.dim,parsPUM.dim);
C = zeros(1,parsPUM.dim);
R = 1;
[ptchPUM, dataXPUM, dataYPUM] = getPtch(C,R,n,parsPUM);
parsPUM.n = n;
parsPUM.M = size(dataYPUM.nodes,1);
parsPUM.N = size(dataXPUM.nodes,1);
%
% Initialise workers
%
delete(gcp('nocreate'));
p = parpool("Processes",8);

%% Serial compute global matrices RBF-PUM
conGlobMat(parsPUM,dataYPUM,dataXPUM,ptchPUM);

%% Parallel compute global matrices RBF-PUM
conGlobMatPar(parsPUM,dataYPUM,dataXPUM,ptchPUM);

%% Serial compute global matrices RBF-FD
conGlobMat(parsFD,dataYFD,dataXFD,ptchFD);

%% Parallel compute global matrices RBF-FD
conGlobMatPar(parsFD,dataYFD,dataXFD,ptchFD);

