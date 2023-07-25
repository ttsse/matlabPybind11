% ----------------------------------------------------------------------------------
% cylPUWeights.m -- compute Partition of Unity weights, their first and
%                   second derivatives. psi = 'wendland_c2' or 'bump'
% Copyright (c) Elisabeth Larsson <eisabeth.larsson@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% ----------------------------------------------------------------------------------

function pu = cylPUWeights(psi,C,R,Z,T,x,plist)

[Np,dim] = size(C);
xloc = globPts2Loc(C,T,x,plist);
%
% Note that gradients are computed in the global coordinate system
%
[phi,phigrad,phiLap,phiH] = cylWeightRBF(psi,C,R,Z,T,xloc,plist);
%
% Compute the sums of the generating functions and their derivatives
%
s = sum(phi,2);

for d=1:dim
  sgrad{d} = sum(phigrad{d},2);
  for k=1:dim
    sH{d,k} = sum(phiH{d,k},2);
  end  
end  

sLap = sum(phiLap,2);
%
% Initialize the data structure as empty and then compute weights
%
[pu(1:Np).w] = deal(zeros(0,1));
for i=1:Np
  loc = find(plist(:,i));
  if (length(loc)>0) % There are points in this patch
    pu(i).w = phi(loc,i)./s(loc);
      
    for d=1:dim
      pu(i).wgrad{d} = (phigrad{d}(loc,i) ...
                        - pu(i).w.*sgrad{d}(loc))./s(loc);
    end
        
    pu(i).wLap = phiLap(loc,i);
    for d=1:dim
      pu(i).wLap = pu(i).wLap -2*pu(i).wgrad{d}.*sgrad{d}(loc);
    end
    pu(i).wLap = (pu(i).wLap - pu(i).w.*sLap(loc))./s(loc);
    %
    % Here we get H and grad matrix product as well, but with the
    % negative sign. Perhaps we need to do it on component form.
    % (Hphi - nabla w^T nabla s + nabla s^T nabla w + wHs)/s 
    %  
    for d=1:dim
      for k=1:dim
        pu(i).wH{d,k} = (phiH{d,k}(loc,i) - ...
                         pu(i).wgrad{d}.*sgrad{k}(loc) - ...
                         pu(i).wgrad{k}.*sgrad{d}(loc) - ...
                         pu(i).w.*sH{d,k}(loc))./s(loc);
      end
    end
  end  
end

  
  
  
  
  