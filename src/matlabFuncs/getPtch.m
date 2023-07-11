%
% Patch generation routine
%
function ptch = getPtch(geom,P,C,R,del)
    dim = size(C,2);
    if strcmp(geom,"ball")
        dimRat = [1 4/pi 8/(4*pi/3)];
        P = ceil(dimRat(dim)*P);
        Rsq = R*(dim)^(1/2);
    elseif strcmp(geom,"cube")
        Rsq = R;
    end
    dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];              % constants for getting cube side length /2
    ptch.R = (2*Rsq)/((2-del)*ceil(P^(1/dim)) - del); % Ensure overlap in cartesian domain (Along diagonal)
    if dim == 1
        Cx = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        ptch.C = Cx(:);
    elseif dim == 2
        Cx = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        Cy = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        [Cx,Cy] = meshgrid(Cx,Cy);
        ptch.C = [Cx(:),Cy(:)];
    elseif dim == 3
        Cx = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        Cy = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        Cz = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        [Cx,Cy,Cz] = meshgrid(Cx,Cy,Cz);
        ptch.C = [Cx(:),Cy(:),Cz(:)];
    end
    ptch.R = ptch.R*ones(size(ptch.C,1),1);
    if strcmp(geom,"ball")
        tol = 1e-6;
        % Remove patches with centres that are outside the domain
        idCentreIn = find(sum((ptch.C - C).^2,2)<=R.^2);
        ptch.C = ptch.C(idCentreIn,:);
        ptch.R = ptch.R(idCentreIn,:);
        % Get closest patch neighbours based on distance
        [idNei,distNei] = rangesearch(ptch.C,ptch.C,(2*ptch.R(1)-ptch.R(1)*del) + ptch.R(1)*tol);
        % Find points of intersection between patches (2 patches in 2D and 4 in 3D)
        edgePts = [];
        midPt = [];
        dirEdgePt = [];
        ptchId = [];
        for i = 1:size(idNei,1)
            neiTemp = idNei{i};
            distTemp = distNei{i};
            % All possible patch combinations in neighbourhood, groups of 2 in 2D and 4 in 3D.
            ptchComb = [nchoosek(neiTemp(2:end),2^(dim-1) - 1)];
            ptchComb = [ones(size(ptchComb,1),1).*neiTemp(1), ptchComb];
            for j = 1:size(ptchComb,1)
                midPtTemp = sum(ptch.C(ptchComb(j,:),:),1)/size(ptchComb,2);
                tempDist = sqrt(sum((ptch.C(ptchComb(j,2:end),:) - ptch.C(ptchComb(j,1),:)).^2,2));
                tempDir = (ptch.C(ptchComb(j,2:end),:) - ptch.C(ptchComb(j,1),:))./tempDist;
                % Computing directions along which the edge point lies (starting from midPt)
                if dim == 2
                    t = (-ptch.C(neiTemp(1),:) + ptch.C(neiTemp(j+1),:))/norm((-ptch.C(neiTemp(1),:) + ptch.C(neiTemp(j+1),:)),2);
                    t = [t(2), -t(1)];
                    t2 = -t;
                    keep = 1;
                else
                    t = cross(tempDir(1,:),tempDir(2,:));
                    keep = ~dot(t,tempDir(end,:));
                    t2 = cross(tempDir(1,:),-tempDir(2,:));
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
        % extension at least a distance of del*ptch.R from the boundary 
        %
        for j = 1:size(edgePts,1)
            if isempty(find(ptPtchList(j,:))) && sqrt(sum((edgePts(j,:)-C).^2)) <= R + del*ptch.R(ptchId(j,1),:)
                r = R + del*ptch.R(ptchId(j,1)); 
                % Direction along which edge point will move when pacthes
                % are equaly extended
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
end