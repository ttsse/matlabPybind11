function [] = importModules
    clear Eigen scipy np
    
    if pyenv().ExecutionMode == "OutOfProcess"
        terminate(pyenv);
    end
    
    if exist("~/miniconda3/envs/TTSSE/bin/python3")
        pyenv(Version="~/miniconda3/envs/TTSSE/bin/python3");
    elseif exist("~/anaconda3/envs/TTSSE/bin/python3")
        pyenv(Version="~/anaconda3/envs/TTSSE/bin/python3");
    else
        error("Install conda (miniconda or anaconda) and generate environment as instructed before running this code")
    end
    
    pyenv("ExecutionMode","OutOfProcess");
    cpp = py.importlib.import_module('cppFuncs');
    np = py.importlib.import_module('numpy');
    scipy = py.importlib.import_module('scipy.sparse');
end

