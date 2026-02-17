% -------------------------------------------------------------------------
% conGlobMatParExp.m -- Performance experiment for parallelisation of
%                       conGlobMat. Result should show two plots of
%                       wall time to compute global evaluation and Laplace
%                       matrices for both the unfitted RBF-PU and RBF-FD
%                       methods in serial and parallel.
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
clearvars;
close all;
%
% Initialise workers
%
delete(gcp('nocreate'));
p = parpool("Processes",8);
%
% Parameter set
%
N = [100 500 1000 1500 2000 2500 3000 3500 4000];
P = [10 50 100 150 200 250 300 350 400];
for i = 1:length(N)
    %
    % Setup RBF-FD parameters
    %
    parsFD = setPars("geom",'ball',"N",N(i),"display",0,"debug",0,"rbfdeg",3,"method","FD");
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
    parsPUM = setPars("geom",'ball',"P",P(i),"display",0,"debug",0,"rbfdeg",3,"method",'PUM');
    n = nchoosek(parsPUM.rbfdeg+parsPUM.dim,parsPUM.dim);
    C = zeros(1,parsPUM.dim);
    R = 1;
    [ptchPUM, dataXPUM, dataYPUM] = getPtch(C,R,n,parsPUM);
    parsPUM.n = n;
    parsPUM.M = size(dataYPUM.nodes,1);
    parsPUM.N = size(dataXPUM.nodes,1);
    %
    % Save number of centre points (degrees of freedom, DOF)
    %
    resPUM(i,1) = length(dataXPUM.nodes);
    resParPUM(i,1) = length(dataXPUM.nodes);
    resFD(i,1) = length(dataXFD.nodes);
    resParFD(i,1) = length(dataXFD.nodes);
    %
    % Save times to compute serial and parallel code
    %
    resPUM(i,2) = timeit(@() conGlobMat(parsPUM,dataYPUM,dataXPUM,ptchPUM));
    resParPUM(i,2) = timeit(@() conGlobMatPar(parsPUM,dataYPUM,dataXPUM,ptchPUM));
    resFD(i,2) = timeit(@() conGlobMat(parsFD,dataYFD,dataXFD,ptchFD));
    resParFD(i,2) = timeit(@() conGlobMatPar(parsFD,dataYFD,dataXFD,ptchFD));
end
%
% Plot times against DOF
%
figure()
plot(resFD(:,1), resFD(:,2), '-o', 'DisplayName', 'Serial');
hold on;
plot(resParFD(:,1), resParFD(:,2), '-o', 'DisplayName', 'Parallel (8 processes)');
xlabel('N');
ylabel('Wall time (seconds)');
title("conGlobMat: RBF-FD matrices")
legend show;
grid on;
figure()
plot(resPUM(:,1), resPUM(:,2), '-o', 'DisplayName', 'Serial');
hold on;
plot(resParPUM(:,1), resParPUM(:,2), '-o', 'DisplayName', 'Parallel (8 processes)');
xlabel('N');
ylabel('Wall time (seconds)');
title("conGlobMat: RBF-PUM matrices")
legend show;
grid on;