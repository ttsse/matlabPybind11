function [tOut] = cppFuncsTest(N,q,display)
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

M = N*q;
%
% Running python "out of process" variables passed are limited to 2GB  
%
tol = .4;
memSize = (M*N*8)*1e-9 + (M*8)*1e-9;

if memSize >= 2 - tol*2 && pyenv().ExecutionMode == "OutOfProcess"
    error("cppFuncTest: System is too large to solve");
end
%
% Solve dense square system
%
A = rand(N,N);
b = rand(N,1);

tic
% xDenseCpp = double(py.cppFuncs.denseSolve(py.numpy.array(A,copy=false,order='F'),py.numpy.array(b,copy=false,order='F')))';
xDenseCpp = double(cpp.denseSolve(np.array(A),np.array(b)));
tCpp = toc;

tic
xDenseMat = A\b;
tMatlab = toc;

disp(['cppFuncTest: Eigen solve time for dense system: ', num2str(tCpp), 's'])
disp(['cppFuncTest: Matlab solve time for dense system: ', num2str(tMatlab), 's'])

if ~setdiff(xDenseMat,xDenseCpp)
    disp("cppFuncTest: Solution incorrect (denseSolve)");
else
    disp("cppFuncTest: Solution OK (denseSolve)");
end
disp(' ');
tOut = [tCpp, tMatlab];

%
% Solve dense LS system
%
A = rand(M,N);
b = rand(M,1);

tic
% xDenseLScpp = double(py.cppFuncs.denseSolve(py.numpy.array(A,copy=false,order='F'),py.numpy.array(b,copy=false,order='F')))';
xDenseLScpp = double(cpp.denseSolve(np.array(A),np.array(b)))';
tCpp = toc;

tic
xDenseLSmat = A\b;
tMatlab = toc;

disp(['cppFuncTest: Eigen solve time for dense LS system: ', num2str(tCpp), 's'])
disp(['cppFuncTest: Matlab solve time for dense LS system: ', num2str(tMatlab), 's'])

if ~setdiff(xDenseLScpp,xDenseLSmat)
    disp("cppFuncTest: Solution incorrect (denseSolve)");
else
    disp("cppFuncTest: Solution OK (denseSolve)");
end
disp(' ');
tOut = [tOut; tCpp, tMatlab];

%
% Solve sparse LS system
%
A = sprand(M,N,0.01,0.1);
b = rand(size(A,1),1);
tic
% Making row, col and data arrays 1D. !Necessary for scipy!
[BI, BJ, BV] = find(A);

row  = np.array(BI - 1).reshape(int32(-1));
col  = np.array(BJ - 1).reshape(int32(-1));
data = np.array(BV,pyargs('dtype', 'float64')).reshape(int32(-1));
pyB = sp.csc_matrix({data,{row, col}}, {int32(size(A,1)), int32(size(A,2))});
tConv = toc;

tic
xSparseLScpp = double(cpp.sparseSolve(pyB,np.array(b,pyargs('dtype', 'float64'))))';
tCpp = toc;

tic
xSparseLSmat = A\b;
tMatlab = toc;

disp(['cppFuncTest: Conversion time to scipy sparse matrix: ', num2str(tConv), 's']);
disp(['cppFuncTest: Eigen solve time for sparse LS system: ', num2str(tCpp), 's']);
disp(['cppFuncTest: Matlab solve time for sparse LS system: ', num2str(tMatlab), 's']);

if ~setdiff(xSparseLScpp,xSparseLSmat)
    disp("cppFuncTest: Solution incorrect (sparseSolve)");
else
    disp("cppFuncTest: Solution OK (sparseSolve)");
end
tOut = [tOut; tCpp, tMatlab];

end