% ----------------------------------------------------------------------------------
% locPts2Ref.m -- Convert coordinates from local coordinate system to 
%                 reference Patch.
% Copyright (c) Elisabeth Larsson <eisabeth.larsson@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% ----------------------------------------------------------------------------------

function xref = locPts2Ref(H,R,Zk,Rk,x)
  dim = size(x,2);
  %
  % We need to take into account that Z is not symmetric around 0
  %
  xref(:,1:dim-1) = R/Rk*x(:,1:dim-1); 
  Hk = diff(Zk);
  Zmid = 0.5*sum(Zk);
  xref(:,dim) = H/Hk*(x(:,dim)-Zmid);