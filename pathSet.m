% -------------------------------------------------------------------------
% pathSet.m -- Sets paths to source
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
addpath("tests/")
addpath(genpath("src/"))
if isfolder("rbfdiff/")
    addpath("rbfdiff")
    addpath("rbfdiff/rbfdir")
    addpath("rbfdiff/rbfqr1D")
    addpath("rbfdiff/rbfqr2D")
    addpath("rbfdiff/rbfqr3D")
    addpath("rbfdiff/utils")
    addpath("rbfdiff/tests")
else
    error("rbfdiff not found, make sure to download rbfdir repository: git clone https://github.com/elisabethl/rbfdiff.git")
end