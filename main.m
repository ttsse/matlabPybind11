% -------------------------------------------------------------------------
% main.m -- Import modules and set parameters used in RBFsolver.
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------

clearvars;
close all;

importModules;
pathSet;

pars.dim = 2;                            % dim = 1,2 or 3
pars.display = 1;                        % Plot solution
pars.geom = 'ball';                      % ball or cube
pars.mode = 'unfitted';                  % fitted, unfitted or collocation
pars.bcMode = 'weak';                    % strong or weak imposition of boundary conditions (only relevant for fitted)
pars.scaling = 1;                        % Include scaling of the unfitted LS problem
pars.mvCentres = 1;                      % Option to have a Y point on top of all X points inside the domain
pars.q = 3;                              % Oversampling (relevant for unfitted and fitted LS methods)
pars.ep = 0.1;                           % Order of 'phs' basis, shape parameter for other bases
pars.phi = 'rbfqr';                      % Choice of basis 'phs', 'mq', 'gs', 'iq', 'rbfqr'
pars.psi = 'w2';                         % Weight function: wendland_c2 or bump (relevant for RBFPUM)
pars.pdeg = -1;                          % Polynomial extension, not relevant for 'rbfqr'
pars.del = 0.4;                          % Overlap between patches

pars.rbfdeg = 4;
pars.memTol = 0.1; 

pars.P = 250; 
pars.extCoeff = 0.5;
% pars.phi = 'phs';
% pars.ep = .1;
pars.method = 'PUM';
pars.N = 15;
pars.prob = 'Poisson';

[l2Error, h] = RBFsolver(pars);