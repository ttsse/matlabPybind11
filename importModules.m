function [cpp, np, sp] = importModules
    clear Eigen scipy np
    
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
    % Importing python libraries
    %
    pyenv("ExecutionMode","OutOfProcess");
    cpp = py.importlib.import_module('cppFuncs');
    np = py.importlib.import_module('numpy');
    sp = py.importlib.import_module('scipy.sparse');
    %
    % Targeted module "warm-up". Required to measure speed of actual
    % operations effectively
    %
    A = [1,0;0,2];
    warmA = np.array(A);
    warmb = np.array([1,1]);
    [row, col, data] = find(A);
    row  = np.array(row - 1).reshape(int32(-1));
    col  = np.array(col - 1).reshape(int32(-1));
    data = np.array(data,pyargs('dtype', 'float64')).reshape(int32(-1));
    warmAsp = sp.csc_matrix({data,{row, col}}, {int32(size(A,1)), int32(size(A,2))});
    cpp.denseSolve(warmA,warmb);
    cpp.sparseSolve(warmAsp,warmb);
    
end

