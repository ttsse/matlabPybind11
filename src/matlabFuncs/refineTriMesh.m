% ----------------------------------------------------------------------------------
% refinrTriMesh.m -- Split triangles until their size is smaller than triA 
% Copyright (c) 2023 Elisabeth Larsson <elisabeth.larsson@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% ----------------------------------------------------------------------------------
function newdata = refineTriMesh(bdata,triA)
  dim = size(bdata.nodes,2);
  tolLarge = 5;

  i1 = [1 2 3];
  i2 = [2 3 1];
  %
  %
  %
  areaList = triArea(bdata);
  tri = bdata.segment;
  nodes = bdata.nodes;
  %
  % Let inner be denoted by zero and outer by 1.
  %
  type = zeros(size(bdata.nodes,1),1);
  type(bdata.outer) = 1;
  %
  % Exclude very large (incorrect) triangles from the refinement
  %
  pos = find(areaList > triA & areaList <= tolLarge*median(areaList));
  triout = tri;
  tri = tri(pos,:);
  triout(pos,:) = [];
  areaList = areaList(pos);
  n0 = size(bdata.nodes,1);
  while length(pos)>0 % Still triangles to refine
    ntri = size(tri,1);
    %
    % Get the edges for all triangles in the list
    %
    e1 = tri(:,i1); e2 = tri(:,i2);
    edge = [e1(:) e2(:)]; % Order is all triangles first side etc
    %
    % Make a list of unique edges, to avoid adding the same point twice
    %
    [newedge,~,newloc] = unique(sort(edge,2),'rows');
    %
    % Because we are pruning the newpts, we need indirection for the pts
    % Here we create indices for the new edge points. Three sides for each
    % triangle.
    %
    midind = zeros(ntri,3);
    for i=1:3
      midind(:,i) = ntri*(i-1)+(1:ntri); % Linear position of new points
      midind(:,i) = newloc(midind(:,i)); % After removing duplicates
    end  
    midind = midind + n0; % Starting after the old points
    %
    % Four new triangles for each old
    %
    newtri = zeros(0,3); % Three corners
    for i=1:dim
      %
      % Triangles with two new midpoints and and one old corner
      %
      newtri = [newtri; midind(:,i) tri(:,i2(i)) midind(:,i2(i))];
    end
    %
    % Triangles with the three new midpoints
    %
    newtri = [newtri; midind];
    %   
    % Here come the new nodes
    %
    x1 = nodes(newedge(:,1),:); x2 = nodes(newedge(:,2),:);
    newnodes = 0.5*(x1+x2);
    %
    % Update the type for each new node
    %
    newtype = type(newedge(:,1));
    type = [type; newtype];
    %
    % The nodes are updated first, so we can use them 
    %
    nodes = [nodes;newnodes];
    %
    % New areas are needed for the new triangles
    %
    x1 = nodes(newtri(:,1),:);
    x2 = nodes(newtri(:,2),:);
    x3 = nodes(newtri(:,3),:);
    x21 = x2-x1;
    x31 = x3-x1;
    newArea = 0.5*sqrt(sum(cross(x21,x31,2).^2,2));
    %
    % Update all values for the next iteration 
    %
    sz = size(newnodes,1);
    n0 = n0 + sz;
    %
    % The refined triangles are removed and the new are added
    %
    tri = newtri;
    areaList = newArea;

    pos = find(areaList > triA);
    tri = tri(pos,:);
    newtri(pos,:) = [];
    triout = [triout;newtri];
    areaList = areaList(pos);
  end
  newdata.nodes = nodes;
  newdata.segment = triout; % Saving all the ready triangles here
  newdata.outer = find(type);
  newdata.inner = find(~type);
