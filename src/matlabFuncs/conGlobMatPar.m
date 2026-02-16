% -------------------------------------------------------------------------
% conGlobMatPar.m -- Compute global approximation matrices for domain and boundary
%                 operators from problem "pars.prob".
% Inputs         -- pars    -> Parameter structure, see setPars. Should
%                              include number of evaluation points, pars.M
%                              and number of local points pars.n
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
%                   pars.M = size(dataY.nodes,1);
%                   pars.N = size(dataX.nodes,1); pars.n = n;
%                   [L,B, Eglobal, Lglobal] = conGlobMatPar(pars,dataY,dataX,ptch);   
%
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function [L,B, Eglobal, Lglobal] = conGlobMatPar(pars,dataY,dataX,ptch)
    % Check input problem
    if ~strcmp(pars.prob,'Poisson')
        error("conGlobMat: Requested problem is not currently implemented")
    end

    xe = dataY.nodes;
    xc = dataX.nodes;
    if strcmp(pars.method,'FD')
        relXc = unique(ptch); % stencil centres that are used directly
        % Construct cell arrays including row / column entry information
        idX = num2cell(knnsearch(xc,xc(relXc,:),'K',pars.n),2);
        [~,evalStencil] = ismember(ptch,relXc);
        idY = accumarray(evalStencil,(1:numel(evalStencil))',[size(idX,1) 1], @(x){x});
        % Initialise index and value cells usedfor sparse assembly
        iRow = cell(pars.n,1); iCol = cell(pars.n,1); 
        Eval = cell(pars.n,1);  Lval = cell(pars.n,1);
        parfor i = 1:length(relXc)
            xcLoc = xc(idX{i},:); % stencil
            Psi = RBFInterpMat(pars.phi,pars.pdeg,pars.ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));
            E = RBFDiffMat(0,Psi,xe(idY{i},:));
            L = RBFDiffMat(1.5,Psi,xe(idY{i},:));
            
            [cc, rr] = meshgrid(idX{i},idY{i});
            iRow{i} = rr(:); iCol{i} = cc(:); 
            Eval{i} = E(:); Lval{i} = L(:);
        end
        iRow = vertcat(iRow{:}); iCol = vertcat(iCol{:});
        Eval = vertcat(Eval{:}); Lval = vertcat(Lval{:});
        Eglobal = sparse(iRow, iCol, Eval, pars.M, pars.N);
        Lglobal = sparse(iRow, iCol, Lval, pars.M, pars.N);
        L = Lglobal(dataY.inner,:);
        B = Eglobal(dataY.bnd,:);
    else
        dim = size(ptch.C,2);
        P = length(ptch.R);
        [w] = weights(pars.psi,1.5,ptch);
        % Initialise index and value cells usedfor sparse assembly
        iE = cell(P,1); jE = cell(P,1); Eval = cell(P,1);
        iL = cell(P,1); jL = cell(P,1); Lval = cell(P,1);
        parfor i = 1:P
            Psi = RBFInterpMat(pars.phi,pars.pdeg,pars.ep,ptch.xc(i).nodes,ptch.C(i,:),ptch.R(i));
            E = RBFDiffMat(0, Psi, ptch.xe(i).nodes);
            B = RBFDiffMat(1, Psi, ptch.xe(i).nodes);
            L = RBFDiffMat(1.5, Psi, ptch.xe(i).nodes);
            rows = ptch.xe(i).globalId; cols = ptch.xc(i).globalId;
            Elocal = w{i}.f'.*E;
            Llocal = zeros(length(rows), length(cols));
            for d = 1:dim
                Llocal = Llocal + 2.*w{i}.grad{d}'.*B{d};
            end
            Llocal = Llocal + w{i}.L'.*E + w{i}.f'.*L;
            [cc, rr] = meshgrid(cols, rows);
            iE{i} = rr(:); jE{i} = cc(:); Eval{i} = Elocal(:);
            iL{i} = rr(:); jL{i} = cc(:); Lval{i} = Llocal(:);
        end
        iE = vertcat(iE{:}); jE = vertcat(jE{:});
        Eval = vertcat(Eval{:});
        iL = vertcat(iL{:}); jL = vertcat(jL{:});
        Lval = vertcat(Lval{:});
        Eglobal = sparse(iE, jE, Eval, pars.M, pars.N);
        Lglobal = sparse(iL, jL, Lval, pars.M, pars.N);
        L = Lglobal(dataY.inner,:);
        B = Eglobal(dataY.bnd,:);
    end
end

