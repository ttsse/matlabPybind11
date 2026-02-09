% -------------------------------------------------------------------------
% movePts.m      -- Move evaluation point (dataY) to closest centre point (dataX). 
%                   Make sure not to move two evaluation points to the same center point.
%                   All center points should be moved to. Interior and then boundary points
%                   moved to avoid changing boundary. Used in RBFsolver.
% Inputs         -- dataX       -> Structure with mandatory fields:
%                   data.nodes  -> Double Nxd array of centre points.
%                   data.inner  -> Double Ninxd array of interior centre points.
%                   data.bnd    -> Double Nbxd array of boundary centre points.
%                   dataY       -> Structure with mandatory fields:
%                   data.nodes  -> Double Mxd array of centre points.
%                   data.inner  -> Double Minxd array of interior centre points.
%                   data.bnd    -> Double Mbxd array of boundary centre points.
% Outputs        -- dataY       -> Updated evaluation point structure.
% Syntax         -- [dataX] = getPts("ball",100,20,[0,0],1,"unfitted",0.5)
%                   [dataY] = getPts("ball",200,20,[0,0],1,"fitted",0)
%                   [dataY] = movePts(dataX,dataY)
%                   
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function dataY = movePts(dataX,dataY)
    %
    % Only move to interior points
    %
    xcIn = dataX.nodes(dataX.inner,:);
    xeIn = dataY.nodes(dataY.inner,:);
    i = 1;
    idY = [1:size(xeIn,1)]'; 
    idX = [1:size(xcIn,1)]'; 
    %
    % Ensure the same Y point is not moved twice and only consider X
    % points that have not been treated
    %
    iXdone = [];             
    iYdone = [];
    while length(iXdone) ~= size(xcIn,1) % ensure all interior center points are treated
        idYmv = knnsearch(xeIn(idY,:),xcIn(idX,:),'k',1);
        [unIdYmv,iX,~] = unique(idYmv,'first');
        xeIn(idY(unIdYmv),:) = xcIn(idX(iX),:);
        iXdone = [iXdone; idX(iX)];
        iYdone = [iYdone; idY(unIdYmv)];
        i = i + 1;
        idY = setdiff(idY,iYdone);
        idX = setdiff(idX,iXdone);
    end
    dataY.nodes(dataY.inner,:) = xeIn;
    %
    % Move boundary points for fitted method only
    %
    if ~isempty(dataX.bnd)
        xcBnd = dataX.nodes(dataX.bnd,:);
        xeBnd = dataY.nodes(dataY.bnd,:);
        i = 1;
        idY = [1:size(xeBnd,1)]'; % All interior evaluation (Y) points
        idX = [1:size(xcBnd,1)]'; % All interior center (X) points
        %
        % Ensure the same Y point is not moved twice and obly consider X
        % points that have not been treated
        %
        iXdone = [];             
        iYdone = [];
        while length(iXdone) ~= size(xcBnd,1) % ensure all interior center points are treated
            idYmv = knnsearch(xeBnd(idY,:),xcBnd(idX,:),'k',1);
            [unIdYmv,iX,~] = unique(idYmv,'first');
            xeBnd(idY(unIdYmv),:) = xcBnd(idX(iX),:);
            iXdone = [iXdone; idX(iX)];
            iYdone = [iYdone; idY(unIdYmv)];
            i = i + 1;
            idY = setdiff(idY,iYdone);
            idX = setdiff(idX,iXdone);
        end
        dataY.nodes(dataY.bnd,:) = xeBnd;
    end
end