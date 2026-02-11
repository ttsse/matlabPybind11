% -------------------------------------------------------------------------
% weights.m      -- Compute partition of unity weights (w) and their derivatives
%                   up to and inluding ndiff order. Note: Uses functions from the 
%                   rbfdiff repository (https://github.com/elisabethl/rbfdiff).
% Inputs         -- psi         -> String. Generating function. w2 for
%                                  Wendland C2 or "bmp" for bump function.
%                   ndiff       -> Positive int. Maximum order of weight
%                                  derivatives to include in w. 0 - projection, 1 - first order
%                                  derivatives, 2- second order derivatives, 1.5 - Laplacian 
%                                  weights to compute.
%                   ptch        -> Patch structure (see getPtch) with fields:
%                   ptch.C      -> Double Pxd array of generated patch centres.
%                   ptch.R      -> Double Px1 array of patch radii. 
%                   ptch.xc             -> Structure with 2xP fields:
%                   ptch.xc(i).nodes    -> Double nxd array of generated
%                                          centre points in patch i.
%                   ptch.xc(i).globalId -> Integer nx1 array of generated
%                                          centre point indices for patch i.
%                   ptch.xe             -> Same as ptch.xc for evaluation points.
% Outputs        -- w           -> 1xP struture with potential fields (depending on ndiff):
%                   w(i).f      -> Double 1xne array of weight function values 
%                                  at ne evaluation points.  
%                   w(i).grad   -> 1xd gradient cell array including double 1xne arrays of 
%                                  derivatives of weight funtion at ne evaluation points.
%                   w(i).L      -> Double 1xne array including values of Laplace 
%                                  of weight function at ne evaluation points. 
%                   w(i).hess   -> dxd Hessian cell array including double 1xne arrays of 
%                                  second derivatives of weight funtion at ne evaluation points.
% Syntax         -- pars = setPars;
%                   n = nchoosek(pars.rbfdeg+pars.dim,pars.dim);
%                   [ptch, dataX, dataY] = getPtch([0,0],1,n,pars);
%                   [w] = weights('bmp',1.5,ptch)
%
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function [w] = weights(psi,ndiff,ptch)
%
% Prepare cell arrays for sums of generating functions and derivatives
%
dim = size(ptch.C,2);
s = zeros(1,length(unique(vertcat(ptch.xe(:).globalId))));
if (ndiff>=1)
    for d = 1:dim
        sGrad{d} = zeros(size(s));
    end
end
if (ndiff==1.5)
    sL = zeros(1,length(unique(vertcat(ptch.xe(:).globalId))));
end
if (ndiff>=2)
    for d = 1:dim
        for k = 1:dim
            sHess{d,k} = zeros(size(s));
        end
    end 
end
%
% Reformulate ndiff to viable operator input for RBFmat function
%
ndiffTot = 0;
op = [];
opDim = [];
ndiffVec = [];
while ndiffTot<ndiff
    [opTemp,opDimTemp,~]=getOpDirect(ndiffTot,dim);
    op = [op,opTemp];                               
    opDim = [opDim, opDimTemp];
    ndiffVec = [ndiffVec, ndiffTot*ones(1,length(opTemp))];
    ndiffTot = ndiffTot + 1;
end
[opTemp,opDimTemp,~]=getOpDirect(ndiff,dim);
op = [op,opTemp];
opDim = [opDim, opDimTemp];
ndiffVec = [ndiffVec, ndiffTot*ones(1,length(opTemp))];
%
% Computegenerating functions, f (bmp or w2) and their requested
% derivatives on each patch.
%
for p = 1:length(ptch.R)
    xLoc = ptch.xe(p).nodes;
    xLoc = (xLoc - ptch.C(p,:))./ptch.R(p); 

    r = xcdist(zeros(1,dim),xLoc,1);    % Euclidean distance matrix
    fTemp = repmat({zeros(size(r))},1,length(op));
    for k=1:length(op)
        fTemp{k} = ((1/ptch.R(p)).^ceil(ndiffVec(k))).*RBFmat(psi,1,r,op{k},opDim{k});
    end
    s(ptch.xe(p).globalId) = s(ptch.xe(p).globalId) + fTemp{1}; % Sum of generating function on all patches
    %
    % Organize the output
    %
    if (ndiff==0)
        f{p}.f = fTemp{1};
    elseif (ndiff==1)
        f{p}.f = fTemp{1};
        f{p}.grad = {};
        for d = 1:dim
            f{p}.grad{end+1} = fTemp{d+1};
            sGrad{d}(ptch.xe(p).globalId) = sGrad{d}(ptch.xe(p).globalId) + fTemp{d+1};
        end
    elseif (ndiff==1.5)
        f{p}.f = fTemp{1};
        f{p}.grad = {};
        for d = 1:dim
            f{p}.grad{end+1} = fTemp{d+1};
            sGrad{d}(ptch.xe(p).globalId) = sGrad{d}(ptch.xe(p).globalId) + fTemp{d+1};
        end
        f{p}.L = fTemp{end};
        sL(ptch.xe(p).globalId) = sL(ptch.xe(p).globalId) + fTemp{end};
    elseif (ndiff==2)
        f{p}.f = fTemp{1};
        f{p}.grad = {};
        for d = 1:dim
            f{p}.grad{end+1} = fTemp{d+1};
            sGrad{d}(ptch.xe(p).globalId) = sGrad{d}(ptch.xe(p).globalId) + fTemp{d+1};
        end
        for k = dim+2:length(op)
            if length(opDim{k}) == 1 
                f{p}.hess{opDim{k},opDim{k}} = fTemp{k};
                sHess{opDim{k},opDim{k}}(ptch.xe(p).globalId) = sHess{opDim{k},opDim{k}}(ptch.xe(p).globalId)+fTemp{k};
            else
                f{p}.hess{opDim{k}(1),opDim{k}(2)} = fTemp{k};
                f{p}.hess{opDim{k}(2),opDim{k}(1)} = fTemp{k};
                sHess{opDim{k}(1),opDim{k}(2)}(ptch.xe(p).globalId) = sHess{opDim{k}(1),opDim{k}(2)}(ptch.xe(p).globalId)+fTemp{k};
                sHess{opDim{k}(2),opDim{k}(1)}(ptch.xe(p).globalId) = sHess{opDim{k}(2),opDim{k}(1)}(ptch.xe(p).globalId)+fTemp{k};
            end
        end
    end
end
%
% Now compute weights using Shepard's method, w = f/s
%
for p = 1:length(ptch.R)
    sLoc = s(ptch.xe(p).globalId);
    if (ndiff==0)
        w{p}.f = f{p}.f./sLoc;
    elseif (ndiff==1)
        w{p}.f = f{p}.f./sLoc;
        w{p}.grad = {};
        for d = 1:dim
            w{p}.grad{end+1} = (f{p}.grad{d} - w{p}.f.*sGrad{d}(ptch.xe(p).globalId))./sLoc;
        end
    elseif (ndiff==1.5)
        w{p}.f = f{p}.f./sLoc;
        w{p}.L = (f{p}.L-w{p}.f.*sL(ptch.xe(p).globalId))./sLoc;
        w{p}.grad = {};
        for d = 1:dim
            w{p}.grad{end+1} = (f{p}.grad{d} - w{p}.f.*sGrad{d}(ptch.xe(p).globalId))./sLoc;
            w{p}.L = w{p}.L - 2.*w{p}.grad{end}.*sGrad{d}(ptch.xe(p).globalId)./sLoc;
        end
    elseif (ndiff==2)
        w{p}.f = f{p}.f./sLoc;
        w{p}.grad = {};
        for d = 1:dim
            w{p}.grad{end+1} = (f{p}.grad{d} - w{p}.f.*sGrad{d}(ptch.xe(p).globalId))./sLoc;
        end
        for k = dim+2:length(op)
            if length(opDim{k}) == 1 
                w{p}.hess{opDim{k},opDim{k}} = (f{p}.hess{opDim{k},opDim{k}}-...
                                               2.*w{p}.grad{opDim{k}}.*sGrad{opDim{k}}(ptch.xe(p).globalId)-...
                                               sHess{opDim{k},opDim{k}}(ptch.xe(p).globalId))./sLoc;
            else
                w{p}.hess{opDim{k}(1),opDim{k}(2)} = (f{p}.hess{opDim{k}(1),opDim{k}(2)}-...
                                                       w{p}.grad{opDim{k}(1)}.*sGrad{opDim{k}(2)}(ptch.xe(p).globalId)-...
                                                       w{p}.grad{opDim{k}(2)}.*sGrad{opDim{k}(1)}(ptch.xe(p).globalId)- ...
                                                       sHess{opDim{k}(1),opDim{k}(2)}(ptch.xe(p).globalId))./sLoc;
                w{p}.hess{opDim{k}(2),opDim{k}(1)} = w{p}.hess{opDim{k}(1),opDim{k}(2)};
            end
        end
    end
end

function [op,opDim,derVec]=getOpDirect(ndiff,dim)
I = eye(dim);
if (ndiff==0)
    op{1} = '0';
    opDim{1} = 1;
    derVec{1} = zeros(1,dim);
elseif (ndiff==1)
    for k=1:dim
        op{k} = '1';
        opDim{k} = k;
        derVec{k} = I(k,:);
    end
elseif (ndiff==1.5)
    op{1} = 'L';
    opDim{1} = dim;
    derVec{1} = 2*I;
elseif (ndiff==2)
    pos = 1;
    for k=1:dim
        op{pos} = '2';
        opDim{pos} = k;
        derVec{pos} = 2*I(k,:);
        pos = pos + 1;
        for j=k+1:dim
            op{pos} = 'm2';
            opDim{pos} = [j k];
            derVec{pos} = I(k,:) + I(j,:);
            pos = pos + 1;
        end    
    end
elseif (ndiff > 2)
    if mod(ndiff,2)==0
        error("weights: Weight derivatives of order higher than 2 are not implemented!")
    end    
end

