% -------------------------------------------------------------------------
% conGlobMat.m -- Compute global approximation matrices for domain and boundary
%                 operators from problem "pars.prob".
% Inputs         -- pars    -> Parameter structure, see setPars.
%                   dataY   -> Centre point structure. See getPts.
%                   dataX   -> Evaluation point structure. See getPts.
%                   ptch    -> if pars.method = 'PUM': patch structure, see
%                              getPtch.
%                           -> if pars.method = 'FD': Double Mx1 array of
%                              closest centre point indices to each  evaluation point.
% Outputs        -- L       -> Double MinxN array approximating Laplacian on interior
%                              domain points, xe.
%                   B       -> Double MbndxN array approximating evaluation
%                              (Dirichlet) operator on boundary points.
%                   Lglobal -> Double MxN array approximating Laplacian on all 
%                              evaluation points, xe.
%                   Eglobal -> Double MxN array approximating evaluation operator
%                              on all evaluation points, xe.
% Syntax         -- pars = setPars;
%                   n = nchoosek(pars.rbfdeg+pars.dim,pars.dim);
%                   [ptch, dataX, dataY] = getPtch([0,0],1,n,pars);
%                   [L,B, Eglobal, Lglobal] = conGlobMat(pars,dataY,dataX,ptch)   
%
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function [L,B, Eglobal, Lglobal] = conGlobMat(pars,dataY,dataX,ptch)
    % Check input problem
    if ~strcmp(pars.prob,'Poisson')
        error("conGlobMat: Requested problem is not currently implemented")
    end

    xe = dataY.nodes;
    xc = dataX.nodes;
    if strcmp(pars.method,'FD')
        relXc = unique(ptch); % stencil centres that are used directly
        Eglobal = spalloc(pars.M,pars.N,pars.M*pars.n);
        Lglobal = spalloc(pars.M,pars.N,pars.M*pars.n);
        for i = 1:length(relXc)
            [idX,~] = knnsearch(xc,xc(relXc(i),:),'K',pars.n);
            xcLoc = xc(idX,:); % stencil
            Psi = RBFInterpMat(pars.phi,pars.pdeg,pars.ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));
            idY = find(ptch==relXc(i));
            E = RBFDiffMat(0,Psi,xe(idY,:));
            L = RBFDiffMat(1.5,Psi,xe(idY,:));
            
            Eglobal(idY,idX) = E + Eglobal(idY,idX);
            Lglobal(idY,idX) = L + Lglobal(idY,idX);
        end
        L = Lglobal(dataY.inner,:);
        B = Eglobal(dataY.bnd,:);
    else
        dim = size(ptch.C,2);
        P = length(ptch.R);
        Eglobal = spalloc(pars.M,pars.N,pars.M*pars.n);
        Lglobal = spalloc(pars.M,pars.N,pars.M*pars.n);
        [w] = weights(pars.psi,1.5,ptch);
        for i = 1:P
            Psi = RBFInterpMat(pars.phi,pars.pdeg,pars.ep,ptch.xc(i).nodes,ptch.C(i,:),ptch.R(i));
            E = RBFDiffMat(0,Psi,ptch.xe(i).nodes);
            B = RBFDiffMat(1,Psi,ptch.xe(i).nodes);
            L = RBFDiffMat(1.5,Psi,ptch.xe(i).nodes);
            Eglobal(ptch.xe(i).globalId,ptch.xc(i).globalId) = w{i}.f'.*E + Eglobal(ptch.xe(i).globalId,ptch.xc(i).globalId);
            for d = 1:dim
                Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId) =  2.*w{i}.grad{d}'.*B{d} + Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId);
            end
            Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId) = w{i}.L'.*E + w{i}.f'.*L + Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId);
        end
        L = Lglobal(dataY.inner,:);
        B = Eglobal(dataY.bnd,:);
    end
end

