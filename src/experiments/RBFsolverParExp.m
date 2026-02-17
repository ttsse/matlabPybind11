% -------------------------------------------------------------------------
% RBFsolverParExp.m --  Performance experiments for entire RBFsolver
%                       routine.
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
clearvars;
close all;
%
% Set modules and parameters
%
[cpp, np, sp] = importModules;
delete(gcp('nocreate'));
parpool("Processes",8);
%
% Parameter set
%
% N = [100 500 1000 1500 2000 2500 3000 3500 4000];
% P = [10 50 100 150 200 250 300 350 400];

N = [200 500 1500 4500 18000];
P = [4 20 120 450 1800];
for i = 1:length(N)
    %
    % Setup parameters
    %
    parsFD = setPars("parOn",0,"geom",'ball',"N",N(i),"display",0,"debug",0,"rbfdeg",4,"method",'FD');   
    parsPUM = setPars("parOn",0,"geom",'ball',"P",P(i),"display",0,"debug",0,"rbfdeg",4,"method",'PUM');
    parsFDpar = setPars("parOn",1,"geom",'ball',"N",N(i),"display",0,"debug",0,"rbfdeg",4,"method",'FD');   
    parsPUMpar = setPars("parOn",1,"geom",'ball',"P",P(i),"display",0,"debug",0,"rbfdeg",4,"method",'PUM');
    n = nchoosek(parsPUMpar.rbfdeg+parsPUMpar.dim,parsPUMpar.dim);
    %
    % Save times to compute serial and parallel code
    %
    tRes.PUM{1}(i) = timeit(@() RBFsolver(parsPUM, cpp, np, sp));
    tRes.PUM{2}(i) = timeit(@() RBFsolver(parsPUMpar, cpp, np, sp));
    tRes.FD{1}(i) = timeit(@() RBFsolver(parsFD, cpp, np, sp));
    tRes.FD{2}(i) = timeit(@() RBFsolver(parsFDpar, cpp, np, sp));
end
%
% Plot times against requested DOF
%
figure()
plot(N, tRes.FD{1}, '-o', 'DisplayName', 'Serial');
hold on;
plot(N,tRes.FD{2}, '-o', 'DisplayName', 'Parallel (8 processes)');
xlabel('N');
ylabel('Wall time (seconds)');
title("RBFsolver: RBF-FD method")
legend show;
grid on;
figure()
plot(P.*n, tRes.PUM{1}, '-o', 'DisplayName', 'Serial');
hold on;
plot(P.*n, tRes.PUM{2}, '-o', 'DisplayName', 'Parallel (8 processes)');
xlabel('N');
ylabel('Wall time (seconds)');
title("RBFsolver: RBF-PU method")
legend show;
grid on;