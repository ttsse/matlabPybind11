function [] = testMemory(N,M,tol)

memSize = N*8*1e-9 + M*8*1e-9;

if memSize >= 2 - tol*2 && pyenv().ExecutionMode == "OutOfProcess"
    error("cppFuncTest: System is too large to solve");
end

end