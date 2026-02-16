% -------------------------------------------------------------------------
% main.m -- Import modules and set parameters, run RBFsolver.
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
clearvars;
close all;
%
% Set path, modules and parameters
%
pathSet;
[cpp, np, sp] = importModules;
pars = setPars;
results = RBFsolver(pars, cpp, np, sp);