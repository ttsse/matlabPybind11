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
parpool;
%
% Parameter set
%
N = [1000 3000 5000 7000 9000 11000];
P = [80 190 300 410 520 630];
for i = 1:length(N)
    %
    % Setup parameters
    %
    parsFD = setPars("cppOn",1,"parOn",0,"geom",'ball',"N",N(i),"display",0,"debug",0,"rbfdeg",4,"method",'FD');   
    parsPUM = setPars("cppOn",1,"parOn",0,"geom",'ball',"P",P(i),"display",0,"debug",0,"rbfdeg",4,"method",'PUM');
    parsFDpar = setPars("cppOn",1,"parOn",1,"geom",'ball',"N",N(i),"display",0,"debug",0,"rbfdeg",4,"method",'FD');   
    parsPUMpar = setPars("cppOn",1,"parOn",1,"geom",'ball',"P",P(i),"display",0,"debug",0,"rbfdeg",4,"method",'PUM');

    %
    % Get number of centre points by running once
    %
    resPUM = RBFsolver(parsPUMpar, cpp, np, sp);
    resFD = RBFsolver(parsFDpar, cpp, np, sp);
    N_PUM(i) = length(resPUM.u);
    N_FD(i) = length(resFD.u);
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
plot(N_FD, tRes.FD{1}, 'b-o', 'DisplayName', 'RBF-FD: Serial', 'LineWidth', 3,"MarkerSize",12);
hold on;
plot(N_FD,tRes.FD{2}, 'r-o', 'DisplayName', 'RBF-FD: Parallel (8 processes)', 'LineWidth', 3,"MarkerSize",12);
plot(N_PUM, tRes.PUM{1}, 'k--o', 'DisplayName', 'RBF-PUM: Serial', 'LineWidth', 3,"MarkerSize",12);
hold on;
plot(N_PUM, tRes.PUM{2}, 'g--o', 'DisplayName', 'RBF-PUM: Parallel (8 processes)', 'LineWidth', 3,"MarkerSize",12);
g = gca;
g.FontSize = 24;
xlim([0,12000]);
ylim([0, 25])
xlabel('N','FontSize',32);
ylabel('wall time (seconds)','FontSize',32);
l = legend;
l.Location = "northwest";
l.FontSize = 24;
grid on;