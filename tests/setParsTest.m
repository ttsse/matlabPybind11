function tests = setParsTest
tests = functiontests(localfunctions);
end
%
% Spurious input
%
function testFunctionOne(testCase)
    input = {"dim",4};
    verifyError(testCase,@() callWithOutput(input),"setPars: ""dim"" should be of the expected data type.");
end

function callWithOutput(input)
    [~] = setPars(input{:});
end

