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
p = parpool;
%
% Parameter set
%
N = [1000 3000 5000 7000 9000 11000];
P = [80 190 300 410 520 630];
for i = 1:length(N)
    %
    % Setup RBF-FD parameters
    %
    parsFD = setPars("geom",'ball',"N",N(i),"display",0,"debug",0,"rbfdeg",4,"method","FD");
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
    parsPUM = setPars("geom",'ball',"P",P(i),"display",0,"debug",0,"rbfdeg",4,"method",'PUM');
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
plot(resFD(:,1), resFD(:,2), 'r-o', 'LineWidth', 3,'DisplayName', 'RBF-FD: Serial',"MarkerSize",12);
hold on;
plot(resParFD(:,1), resParFD(:,2), 'b-o', 'LineWidth', 3, 'DisplayName', 'RBF-FD: Parallel (8 processes)',"MarkerSize",12);
plot(resPUM(:,1), resPUM(:,2), 'k--o', 'LineWidth', 3, 'DisplayName', 'RBF-PUM: Serial',"MarkerSize",12);
hold on;
plot(resParPUM(:,1), resParPUM(:,2), 'g--o', 'LineWidth', 3, 'DisplayName', 'RBF-PUM: Parallel (8 processes)',"MarkerSize",12);
g = gca;
g.FontSize = 24;
xlim([0,11000]);
ylim([0, 25])
xlabel('N','FontSize',32);
ylabel('wall time (seconds)','FontSize',32);
l = legend;
l.Location = "northwest";
l.FontSize = 24;
grid on;