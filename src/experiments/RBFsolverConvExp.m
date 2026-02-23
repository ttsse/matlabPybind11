% -------------------------------------------------------------------------
% RBFsolverConvExp.m --  Convergence results for RBFsolver
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
    parsFD = setPars("cppOn",1,"parOn",1,"geom",'ball',"N",N(i),"display",0,"debug",0, ...
                     "rbfdeg",4,"method",'FD',"mode",'unfitted');   
    parsPUM = setPars("cppOn",1,"parOn",1,"geom",'ball',"P",P(i),"display",0,"debug",0,...
                      "rbfdeg",4,"method",'PUM',"mode",'unfitted');
    n = nchoosek(parsPUM.rbfdeg+parsPUM.dim,parsPUM.dim);
    %
    % Save results
    %
    Res.PUM{i} = RBFsolver(parsPUM, cpp, np, sp);
    Res.FD{i} = RBFsolver(parsFD, cpp, np, sp);
end
%
% Collect fill distance and error
%
tmp = [Res.PUM{:}];
errorPUM = [tmp.l2Error];
tmp = [Res.FD{:}];
errorFD = [tmp.l2Error];
tmp = [Res.PUM{:}];
N_PUM = arrayfun(@(s) length(s.u),tmp);
tmp = [Res.FD{:}];
N_FD = arrayfun(@(s) length(s.u),tmp);
%
% Plot convergence for both methods
%
figure();
pFD = polyfit(log(N_FD),log(errorFD),1);
pPUM = polyfit(log(N_PUM),log(errorPUM),1);
nameFD =  ['RBF-FD',',    slope = ',num2str(pFD(1))];
namePUM = ['RBF-PUM',   ', slope = ',num2str(pPUM(1))];
loglog(N_FD,errorFD,'r-o',"LineWidth",3,'DisplayName',nameFD,"MarkerSize",12); hold on
loglog(N_PUM,errorPUM,'b-o',"LineWidth",3,'DisplayName',namePUM,"MarkerSize",12); 
fig = gca;
fig.FontSize = 24;
fig.XMinorGrid = "off";
fig.YMinorGrid = "off";
fig.XLim = [8e2 2e4];
fig.YLim = [1e-6 3e-3];
xlabel('N','FontSize',32);
ylabel('error','FontSize',32);
l = legend;
l.Location = "southwest";
l.FontSize = 24;
plot(N_FD,(N_FD.^pFD(1)).*exp(pFD(2)),'r--','HandleVisibility','off','LineWidth',2)
plot(N_PUM,(N_PUM.^pPUM(1)).*exp(pPUM(2)),'b--','HandleVisibility','off','LineWidth',2)
grid on
grid minor

