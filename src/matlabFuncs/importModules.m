% -------------------------------------------------------------------------
% importModules.m -- Initialise python environment in matlab 
%                    (should correspond to python in conda env). Then
%                    import relevant libraries (including ones made from
%                    pybind). This function can be adjusted to include
%                    more relevant (available) libraries to any implementation.
%                    The function also initializes all libraries to remove
%                    startup overhead.
% Syntax          -- [cpp, np, sp] = importModules;
%
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function [varargout] = importModules
    %
    % Start a new process
    %
    if pyenv().ExecutionMode == "OutOfProcess"
        terminate(pyenv);
    end
    %
    % Import python version from conda environment
    %
    [~,cmdOut] = system('conda info --base');
    cmdOut = strtrim(cmdOut);
    if strcmp(computer,"PCWIN64") % GLNXA64 - Linux, MACI64/MACA64 - macOS / silicon
        pyPath = '\envs\TTSSE\python.exe';
    else
        pyPath = '/envs/TTSSE/bin/python3';
    end
    if exist([cmdOut, pyPath])
        pyenv(Version=[cmdOut, pyPath]);
    else
        error("cppFuncTest: Install conda (miniconda or anaconda) and generate environment as instructed before running this code")
    end
    %
    % Import python libraries
    %
    pyenv("ExecutionMode","OutOfProcess");
    varargout{1} = py.importlib.import_module('cppFuncs');
    varargout{2} = py.importlib.import_module('numpy');
    varargout{3} = py.importlib.import_module('scipy.sparse');
    %
    % Targeted module "warm-up". Required to measure speed correctly.
    %
    A = [1,0;0,2];
    warmA = varargout{2}.array(A);
    warmb = varargout{2}.array([1,1]);
    [row, col, data] = find(A);
    row  = varargout{2}.array(row - 1).reshape(int32(-1));
    col  = varargout{2}.array(col - 1).reshape(int32(-1));
    data = varargout{2}.array(data,pyargs('dtype', 'float64')).reshape(int32(-1));
    warmAsp = varargout{3}.csc_matrix({data,{row, col}}, {int32(size(A,1)), int32(size(A,2))});
    varargout{1}.denseSolve(warmA,warmb);
    varargout{1}.sparseSolve(warmAsp,warmb);
end

