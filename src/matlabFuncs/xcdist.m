% ----------------------------------------------------------------------------------
% xcdist.m -- Compute distance between evaluation and collocation points,
%             r=||xi-cj|| together with componentwise signed differences.
% Copyright (c) Elisabeth Larsson <elisabeth.larsson@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% ----------------------------------------------------------------------------------

function r=xcdist(x,c,all)

if (nargin==1)
  all=0;
  c=x;
elseif (nargin==2)
  all=0;
end
  
[np,nd] = size(x);
nc = size(c,1);
nr = 1 + all*nd;

% r contains r and each ri, i=1...nd
r = zeros(np,nc,nr);
if (np==0 | nc==0)
  return
end  
for d=1:nd
  [pi,pj] = meshgrid(c(:,d),x(:,d));
  r(:,:,1) = r(:,:,1) + (pi-pj).^2;
  if (all)
    r(:,:,d+1) = pj-pi;
  end
end
r(:,:,1) = sqrt(r(:,:,1));






