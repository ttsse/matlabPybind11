% ----------------------------------------------------------------------------------
% globPts2Loc.m -- Convert coordinates from global coordinate system to local.
% Copyright (c) Elisabeth Larsson <eisabeth.larsson@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% ----------------------------------------------------------------------------------

function xloc = globPts2Loc(C,T,xglob,plist)
  for k=1:length(T)
    if (nargin==4)  
      if (isstruct(plist))
        locpts = plist(k).ind;
      else
        locpts = find(plist(:,k));
      end  
      xg = xglob(locpts,:); 
    else % Already a structure
      xg = xglob(k).pts;
    end
    xg = xg - ones(size(xg,1),1)*C(k,:);
    xg = xg*T(k).Q;
    
    xloc(k).pts = xg;
  end
  
    
    
    