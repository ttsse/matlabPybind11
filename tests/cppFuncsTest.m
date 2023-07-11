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
if exist("~/miniconda3/envs/TTSSE/bin/python3")
    pyenv(Version="~/miniconda3/envs/TTSSE/bin/python3");
elseif exist("~/anaconda3/envs/TTSSE/bin/python3")
    pyenv(Version="~/anaconda3/envs/TTSSE/bin/python3");
else
    error("cppFuncTest: Install conda (miniconda or anaconda) and generate environment as instructed before running this code")
end
%
% Importing python libraries
%
pyenv("ExecutionMode","OutOfProcess");
cpp = py.importlib.import_module('cppFuncs');
np = py.importlib.import_module('numpy');
scipy = py.importlib.import_module('scipy.sparse');

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
xDenseCpp = double(py.cppFuncs.denseSolve(py.numpy.array(A,copy=false,order='F'),py.numpy.array(b,copy=false,order='F')))';
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
xDenseLScpp = double(py.cppFuncs.denseSolve(py.numpy.array(A,copy=false,order='F'),py.numpy.array(b,copy=false,order='F')))';
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
A = sprand(M,N,0.01);
tic
[BI, BJ, BV] = find(A);
pyB = py.scipy.sparse.csc_matrix({BV, {uint64(BI-1) uint64(BJ-1)}}, {uint64(size(A,1)), uint64(size(A,2))});
tConv = toc;

tic
xSparseLScpp = double(py.cppFuncs.sparseSolve(pyB,py.numpy.array(b,copy=false,order='F')))';
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