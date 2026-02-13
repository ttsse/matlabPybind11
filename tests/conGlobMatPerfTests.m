% -------------------------------------------------------------------------
% conGlobMatPerfTests.m -- General performance tests for parallelisation of
%                         conGlobMat
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
clearvars;
close all;

pathSet; 
pars = setPars("N",100,"P",100,"display",0,"debug",0,"rbfdeg",5,"method","FD");
% Get centre points (xc)
n = nchoosek(pars.rbfdeg+pars.dim,pars.dim);
C = zeros(1,pars.dim);
if strcmp(pars.geom,'cube')
    R = 1*(pars.dim)^(1/2);
elseif strcmp(pars.geom,'ball')
    R = 1;
end
[dataX] = getPts(pars.geom,pars.N,n,C,R,pars.mode,pars.extCoeff);
xc = dataX.nodes;         
pars.M = ceil(pars.q*pars.N);   % Number of evaluation points (Y)
[dataY] = getPts(pars.geom,pars.M,n,C,R,"fitted",0);
if pars.mvCentres 
    dataY = movePts(dataX,dataY);
end
xe = dataY.nodes;
xc = xc(unique(knnsearch(xc,xe,'K',ceil(pars.extCoeff*n+pars.ep))),:);
pars.N = size(xc,1);
dataX.nodes = xc;
ptch = knnsearch(xc,xe,'K',1); % This is a list of closest centre points to eval points (centres)
pars.n = n;
pars.M = size(dataY.nodes,1);
pars.N = size(dataX.nodes,1); 
% delete(gcp('nocreate'));
% p = parpool("Processes",8);

%% Serial compute global matrices
conGlobMat(pars,dataY,dataX,ptch);

%% Parallel compute global matrices
% conGlobMatPar(pars,dataY,dataX,ptch);

