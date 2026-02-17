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
parpool("Processes",8);
%
% Parameter set
%
N = [200 500 1500 4500 18000];
P = [4 20 120 450 1800];
for i = 1:length(N)
    %
    % Setup parameters
    %
    parsFD = setPars("parOn",1,"geom",'ball',"N",N(i),"display",0,"debug",0, ...
                     "phi",'phs',"pdeg",4,"ep",3,"method",'FD');   
    parsPUM = setPars("parOn",1,"geom",'ball',"P",P(i),"display",0,"debug",0,...
                      "rbfdeg",4,"method",'PUM');
    n = nchoosek(parsPUM.rbfdeg+parsPUM.dim,parsPUM.dim);
    %
    % Save times to compute serial and parallel code
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
hPUM = [tmp.h];
tmp = [Res.FD{:}];
hFD = [tmp.h];
%
% Plot convergence for both methods
%
figure();
c = colormap("turbo");
nameFD = ['unfitted RBF-FD,    order = ',num2str(pFD(1))];
namePUM = ['unfitted RBF-PUM, order = ',num2str(pPUM(1))];
loglog(hFD,errorFD,'r-o',"LineWidth",2,'DisplayName',nameFD); hold on
loglog(hPUM,errorPUM,'b-o',"LineWidth",2,'DisplayName',namePUM); 
fig = gca;
fig.FontSize = 16;
fig.XMinorGrid = "off";
fig.YMinorGrid = "off";
fig.XLim = [0.01 1];
fig.YLim = [1e-7 1];
pFD = polyfit(log(hFD),log(errorFD),1);
pPUM = polyfit(log(hPUM),log(errorPUM),1);
l = legend;
l.Location = "southeast";
l.FontSize = 18;
plot(hFD,(hFD.^pFD(1)).*exp(pFD(2)),'r--','HandleVisibility','off','LineWidth',1.5)
plot(hPUM,(hPUM.^pPUM(1)).*exp(pPUM(2)),'b--','HandleVisibility','off','LineWidth',1.5)
grid on

