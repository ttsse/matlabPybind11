%
% Make sure not to move two evaluation points to the same center point. All
% center points should be moved to. Move interior and then boundary points
% to avoid changing boundary
%
function dataY = movePts(dataX,dataY)
    % Move interior points
    xcIn = dataX.nodes(dataX.inner,:);
    xeIn = dataY.nodes(dataY.inner,:);
    i = 1;
    idY = [1:size(xeIn,1)]'; % All interior evaluation (Y) points
    idX = [1:size(xcIn,1)]'; % All interior center (X) points
    %
    % Ensure the same Y point is not moved twice and obly consider X
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