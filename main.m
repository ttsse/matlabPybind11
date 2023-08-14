% -------------------------------------------------------------------------
% main.m -- simulate derfomation of 3D thin, linearly elastic bodies
%           using the RBF-PU method
% Copyright (c) 2023 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------

clearvars;
close all;
clear Eigen scipy np

if pyenv().ExecutionMode == "OutOfProcess"
    terminate(pyenv);
end

pyenv(Version="~/anaconda3/envs/TTSSE3/bin/python3");
pyenv("ExecutionMode","OutOfProcess");
cpp = py.importlib.import_module('cppFuncs');
np = py.importlib.import_module('numpy');
scipy = py.importlib.import_module('scipy.sparse');

warning on;

addpath(genpath('src/'))
addpath(genpath('tests/'))

mode = 'manufactured'; % manufactured or physical 
man_sol = 'gauss'; % gauss cosine 
% Material Parameters
pars.E = 1e6; % Young's Modulus
pars.nu = 0.3; % Poisson's ratio
pars.rho = 0;
pars.g = 9.81;
p_max = -30; % Boundary condition for cuboid

numPatches = [20];
nLoc = [56]; 
pars.oversamp = 3;
pars.ep = 0.1;%0.5; 
pars.overlap = 0.05;
pars.dataTol = 1; 
pars.Hptch = 0.65;
pars.Rptch = 1;
pars.errmode = 'exact';
pars.np0 = numPatches(1);
flags.visualize = 0;

for i = 1:length(nLoc)
    for j = 1:length(numPatches)
    pars.nLoc = nLoc(i);
    pars.numPatches = numPatches(j);
    l = 200; w = 50; h = 5;
    [bdata, geom] = makeFEMRectangle(l,w,h,0);
    fldr = 'src/experiments/cuboid/';
    
    flags.iteration = j;
    
    matfile = ['Cover_cuboid_Patches_' num2str(pars.numPatches) '_np0_' ...
           num2str(pars.np0) 'nLoc_' num2str(pars.nLoc) 'oversamp_' num2str(pars.oversamp) ...
           'overlap_' num2str(pars.overlap) 'ep_' num2str(pars.ep) '.mat'];
    if exist(['src/data/cuboid/' matfile], 'file')
        load(['src/data/cuboid/' matfile]);
        pars.Hscale = Hscale;
        pars.Rscale = Rscale;
    else
        error("The chosen discretization is not yet available!")
    end

    pars.mode = mode;
    switch mode 
        case 'manufactured'
            syms x1 x2 x3
            switch man_sol     
                case 'cosine'
                    u_exact{1} = [0.*x1 -(x2/w).^3 -h.*cos((pi/l).*x1)];
                    pars.break_points = [];
                case 'gauss'  
                    dip = 5;
                    B_edge = ((l/2)^2+(w/2)^2-1);
                    C_1 = (dip)/(exp(B_edge/1e4)-exp((-1)/1e4));
                    C_2 = -C_1*exp(B_edge/1e4);
                    u_exact{1} = [0.*x1 -(x2/w)^3 C_1*exp((-1)./(1e4)*(1-(x1.*x1)-(x2.*x2))) + C_2];
                    pars.break_points = [];
            end
            pars.u_exact = u_exact;
            pars.man_sol = man_sol;
        case 'physical'
            pars.body_force.x = 0;
            pars.body_force.y = 0;
            pars.body_force.z = -pars.rho*pars.g;  
            pars.pressure_bottom = 0;                          
            pars.pressure_top = p_max;
    end
    
    switch mode 
        case 'manufactured'
            exper = ['3D_cuboid_nLoc' num2str(pars.nLoc) '_np0' ...
            num2str(pars.np0) '_Patches' ...
            num2str(numPatches(j)) ...
            '_delta' num2str(pars.overlap) '_q' ...
            num2str(pars.oversamp) '_ep' num2str(pars.ep) '_nu' num2str(pars.nu*100) '_manufactured_'  man_sol];
        case 'physical'
            exper = ['3D_cuboid_nLoc' num2str(pars.nLoc) '_np0' ...
            num2str(pars.np0) '_Patches' ...
            num2str(pars.numPatches) ...
            '_delta' num2str(pars.overlap) '_q' ...
            num2str(pars.oversamp) '_nu' num2str(pars.nu*100) '_', bc, '_p' ...
            num2str(-p_max) '_rho' num2str(round(pars.rho))];
    end
    
    pars.scaling.dirichlet = @(h,hy,dim) ((hy^((dim-1)/2))/(h^(1.5)));
    pars.scaling.neumann = @(h,hy,dim) ((hy^((dim-1)/2))/(h^(0.5)));
    pars.scaling.interior = @(h,hy,dim) (hy^(dim/2));
                
    [pars, Result] = solveElasticity(geom,X_points,Y_points,Patches,pars);
    
    mkdir(fldr)
    save([fldr exper '.mat'], 'Result', 'pars', 'bdata');
    end
end      