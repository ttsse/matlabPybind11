% -------------------------------------------------------------------------
% getPtch.m      -- Patch (lines, circles, spheres) and point generation routine.
% Inputs         -- C                   -> Double 1xd array of domain centre in 
%                                          d dimensions.
%                   R                   -> Double. Domain radius.
%                   n                   -> Positive int. Number of local centre 
%                                          points to generate.
%                   pars                -> Structure with mandatory parameters
%                                          geom, P, del, mode, q, mvCentres. For
%                                          details see setPars.
% Outputs        -- ptch                -> Structure with fields:
%                   ptch.C              -> Double Pxd array of generated patch centres.
%                   ptch.R              -> Double Px1 array of patch radii. 
%                   ptch.xc             -> Structure with 2xP fields:
%                   ptch.xc(i).nodes    -> Double nxd array of generated
%                                          centre points in patch i.
%                   ptch.xc(i).globalId -> Integer nx1 array of generated
%                                          centre point indices for patch i.
%                   ptch.xe             -> Same as ptch.xc for evaluation points.
%                   dataX               -> Centre point structure with fields:
%                   dataX.nodes         -> Double Nxd array of generated points.
%                   dataX.inner         -> Double Ninxd array of interior points.
%                   dataX.outer         -> Double Noutxd array of outer points (empty for "fitted").
%                   dataX.bnd           -> Double Nbxd array of boundary points.
%                   dataY               -> Same as dataX for evaluation points.
% Syntax         -- pars = setPars;
%                   n = nchoosek(pars.rbfdeg+pars.dim,pars.dim);
%                   [ptch, dataX, dataY] = getPtch([0,0],1,n,pars);
%                   plotPtch(ptch,'ball',[0,0],1);
%
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function [ptch,dataX,dataY] = getPtch(C,R,n,pars)
    dim = size(C,2);
    %
    % Radius of Cartesian grid
    %
    if strcmp(pars.geom,"ball")
        dimRat = [1 4/pi 8/(4*pi/3)];
        P = ceil(dimRat(dim)*pars.P);
        Rsq = R*(dim)^(1/2);
    elseif strcmp(pars.geom,"cube")
        Rsq = R;
        P = pars.P;
    end
    dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];              % constants for getting cube side length /2
    ptch.R = (2*Rsq)/((2-pars.del)*ceil(P^(1/dim)) - pars.del); % Ensure overlap in Cartesian grid (along diagonal)
    %
    % Place centres on grid 
    %
    if dim == 1
        Cx = linspace(-1+dimLCoeff(dim)*(1-pars.del)*ptch.R,1-dimLCoeff(dim)*(1-pars.del)*ptch.R,ceil(P^(1/dim)));
        ptch.C = Cx(:);
    elseif dim == 2
        Cx = linspace(-1+dimLCoeff(dim)*(1-pars.del)*ptch.R,1-dimLCoeff(dim)*(1-pars.del)*ptch.R,ceil(P^(1/dim)));
        Cy = linspace(-1+dimLCoeff(dim)*(1-pars.del)*ptch.R,1-dimLCoeff(dim)*(1-pars.del)*ptch.R,ceil(P^(1/dim)));
        [Cx,Cy] = meshgrid(Cx,Cy);
        ptch.C = [Cx(:),Cy(:)];
    elseif dim == 3
        Cx = linspace(-1+dimLCoeff(dim)*(1-pars.del)*ptch.R,1-dimLCoeff(dim)*(1-pars.del)*ptch.R,ceil(P^(1/dim)));
        Cy = linspace(-1+dimLCoeff(dim)*(1-pars.del)*ptch.R,1-dimLCoeff(dim)*(1-pars.del)*ptch.R,ceil(P^(1/dim)));
        Cz = linspace(-1+dimLCoeff(dim)*(1-pars.del)*ptch.R,1-dimLCoeff(dim)*(1-pars.del)*ptch.R,ceil(P^(1/dim)));
        [Cx,Cy,Cz] = meshgrid(Cx,Cy,Cz);
        ptch.C = [Cx(:),Cy(:),Cz(:)];
    end
    ptch.R = ptch.R*ones(size(ptch.C,1),1);
    %
    % If geometry is a ball remove patches outside the domain
    %
    if strcmp(pars.geom,"ball") && dim >= 2
        tol = 1e-6;
        % Remove patches with centres that are outside the domain
        idCentreIn = find(sum((ptch.C - C).^2,2)<=R.^2);
        ptch.C = ptch.C(idCentreIn,:);
        ptch.R = ptch.R(idCentreIn,:);
        % Get closest patch neighbours based on distance
        [idNei,~] = rangesearch(ptch.C,ptch.C,(2*ptch.R(1)-ptch.R(1)*pars.del) + ptch.R(1)*tol);
        % Find points of intersection between patches (interesction of 2 patches in 2D and 4 in 3D)
        edgePts = [];
        midPt = [];
        ptchId = [];
        for i = 1:size(idNei,1)
            neiTemp = idNei{i};
            % Get all possible patch combinations in neighbourhood, groups of 2 in 2D and 4 in 3D.
            ptchComb = [nchoosek(neiTemp(2:end),2^(dim-1) - 1)];
            ptchComb = [ones(size(ptchComb,1),1).*neiTemp(1), ptchComb];
            for j = 1:size(ptchComb,1)
                midPtTemp = sum(ptch.C(ptchComb(j,:),:),1)/size(ptchComb,2); % mid-point between patch centres in neighbourghood
                tempDist = sqrt(sum((ptch.C(ptchComb(j,2:end),:) - ptch.C(ptchComb(j,1),:)).^2,2));
                tempDir = (ptch.C(ptchComb(j,2:end),:) - ptch.C(ptchComb(j,1),:))./tempDist;
                % Computing directions (t, t2) along which the intersection (edge) point lies (starting from midPt)
                if dim == 2
                    t = (-ptch.C(neiTemp(1),:) + ptch.C(neiTemp(j+1),:))/norm((-ptch.C(neiTemp(1),:) + ptch.C(neiTemp(j+1),:)),2);
                    t = [t(2), -t(1)];
                    t2 = -t;
                    keep = 1;
                else
                    t = cross(tempDir(1,:),tempDir(2,:));
                    t2 = cross(tempDir(1,:),-tempDir(2,:));
                    keep = ~dot(t,tempDir(end,:));
                end
                if keep % In 3D some combinations of patch centres do not form a plane
                    tempL = sqrt((ptch.R(1)^2)-(max(tempDist)/2)^2);
                    edgePts = [edgePts; tempL*t + midPtTemp];
                    edgePts = [edgePts; tempL*t2 + midPtTemp];
                    ptchId = [ptchId; ptchComb(j,:); ptchComb(j,:)];
                    midPt = [midPt; midPtTemp; midPtTemp];
                end                    
            end
        end
        [edgePts ,idKeep]= uniquetol(edgePts,tol,'ByRows',true);
        ptchId = ptchId(idKeep,:);
        midPt = midPt(idKeep,:);
        %
        % Find edge points that are not covered
        %
        for i = 1:size(edgePts,1)
            ptPtchList(i,:) = sqrt(sum((edgePts(i,:) - ptch.C).^2,2)) <= ptch.R - ptch.R.*(tol); % Points covered by which patches
        end
        %
        % Extend patch sizes to overlap these points. As a result the
        % extension is at least a distance of pars.del*ptch.R from the boundary 
        %
        for j = 1:size(edgePts,1)
            if isempty(find(ptPtchList(j,:))) && sqrt(sum((edgePts(j,:)-C).^2)) <= R + pars.del*ptch.R(ptchId(j,1),:)
                r = R + pars.del*ptch.R(ptchId(j,1)); 
                % Direction along which edge point will move when patch sizes increase
                gamma = (-midPt(j,:)+edgePts(j,:))./sqrt(sum((edgePts(j,:)-midPt(j,:)).^2)); 
                theta = acos(dot(((C-edgePts(j,:))./sqrt(sum((edgePts(j,:)-C).^2))),gamma));
                if abs(theta) <= pi + tol && abs(theta) >= pi - tol % In case gamma direction is towards the centre of the domain
                    theta = pi;
                end
                SQ = (cos(theta).*sqrt(sum((edgePts(j,:)-C).^2))).^2 - sum((edgePts(j,:)-C).^2) + r^2;
                dMv = sqrt(sum((edgePts(j,:)-C).^2))*cos(theta) + sqrt(SQ);
                y = edgePts(j,:) + dMv.*gamma;
                ptch.R(ptchId(j,:)) = sqrt(sum((y - ptch.C(ptchId(j,1),:)).^2));
            end
        end
    end
    P = length(ptch.R);
    %
    % Get centres (xc)
    %
    if strcmp(pars.mode,"unfitted")
        dataX.nodes = [];            % Centers
        for i = 1:P
            [ptch.xc(i)] = getPts("ball",n,0,ptch.C(i,:),ptch.R(i),"unfitted",0);
            xCglobalId{i} = [(i-1)*n+1:i*n]';
            dataX.nodes = [dataX.nodes; ptch.xc(i).nodes];
        end
        xc = dataX.nodes;
        [ptch.xc.globalId] = xCglobalId{:};
        if strcmp(pars.geom,"ball")
            dataX.inner = find(sqrt(sum((dataX.nodes-C).^2,2))<=R);
        elseif strcmp(pars.geom,"cube")
            dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];  
            dataX.inner = find(all(abs(dataX.nodes-C)<=dimLCoeff(dim)*R,2));
        end
        dataX.outer = setdiff(1:size(dataX.nodes,1),dataX.inner);
        dataX.bnd = [];
    elseif strcmp(pars.mode,"collocation")
        dataX = getPts(pars.geom,n*P,0,C,R,"fitted",0);
        for i = 1:P
            ptch.xc(i).globalId = find(sqrt(sum((dataX.nodes - ptch.C(i,:)).^2,2)) <= ptch.R(i));
            ptch.xc(i).nodes = dataX.nodes(ptch.xc(i).globalId,:);
        end
    elseif strcmp(pars.mode,"fitted") 
        dataX = getPts(pars.geom,n*P,0,C,R,"fitted",0);
        for i = 1:P
            ptch.xc(i).globalId = find(sqrt(sum((dataX.nodes - ptch.C(i,:)).^2,2)) <= ptch.R(i));
            ptch.xc(i).nodes = dataX.nodes(ptch.xc(i).globalId,:);
        end
    end
    %
    % Get evaluation points (xe)
    %
    q = max(pars.q*double(~strcmp(pars.mode,"collocation")),1);
    M = ceil(n*P*q);  
    if ~strcmp(pars.mode,"collocation")
        dataY = getPts(pars.geom,M,0,C,R,"fitted",0);
        %
        % Move evaluation points inside (and on the boundary) to the closest center point.
        %
        if pars.mvCentres && ~strcmp(pars.mode,"collocation")
            dataY = movePts(dataX,dataY);
        end
        xe = dataY.nodes;
        for i = 1:P
            ptch.xe(i).globalId = find(sqrt(sum((xe - ptch.C(i,:)).^2,2)) <= ptch.R(i));
            ptch.xe(i).nodes = xe(ptch.xe(i).globalId,:);
        end
    else
        dataY = dataX;
        ptch.xe = ptch.xc;
    end
end