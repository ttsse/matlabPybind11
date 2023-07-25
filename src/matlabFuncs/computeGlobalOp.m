% ----------------------------------------------------------------------------------
% computeGlobalOp.m -- compute approximate differentiation operator, op,
%                      in a patch covered domain. 
%                      Options for op: 'H' - Hessian, 'J' - Jacobian, 
%                      'L' - Laplace, '0' - Interpolation  
% Copyright (c) 2021 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% ----------------------------------------------------------------------------------

function Eout=computeGlobalOp(approx,C,R,Z,T,op,x,plist,pu)
%
% First we allocate the global sparse matrix.
%
  dim = size(x,2);
  phi = approx(1).phi;
  ep = approx(1).ep;
  pdeg = approx(1).pdeg;
  xc = approx(1).xc; 
  nc = size(xc,1);
  np = size(plist,2);
  N = size(x,1); 
  npt = sum(plist);
  
  if (op == 'H')
      for d_1 = 1:dim
        for d_2 = 1:dim 
            Eout{d_1,d_2} = spalloc(N,nc*np,nc*sum(npt));
        end
      end 
  elseif (op == 'J')
      for d = 1:dim
        Eout{d} = spalloc(N,nc*np,nc*sum(npt));
      end 
  else
      Eout = spalloc(N,nc*np,nc*sum(npt));
  end

  %
  % Then we loop through the patches to create the contributions.
  %
  xloc = globPts2Loc(C,T,x,plist);
  for k=1:np
    %
    % Take the local points in this patch to the reference patch
    %
    xref = locPts2Ref(approx(1).H,approx(1).R,Z(k,:),R(k),xloc(k).pts);
    
    ind = find(plist(:,k));
    nn = npt(k);
    
    re = xcdist(xref,xc,1);
    if (strcmp(phi,'rbfqr'))
      if (dim==2)
        [B,approx(1).Psi] = RBF_QR_diffmat_2D('1',xref,approx(1).Psi);
      else
        [B,approx(1).Psi] = RBF_QR_diffmat_3D('0',xref,approx(1).Psi);
      end
    else
      B = RBFmat(phi,ep,re,'0');
    end  
    %B = RBFMother(phi,ep,re,'0', dim, xref, xc);
    
    if (pdeg >= 0)
      P = polyMat(xref,pdeg);
      B = [B P];
    end
    E = (B/approx(1).U)/approx(1).L;
    E(:,approx(1).piv) = E;
    E = E(:,1:nc);
    
    if (op~='0')
      %
      % Gradient, and make sure to transform matrices to the global
      % coordinate system.
      %
      dscale = [approx(1).R/R(k)*ones(1:dim-1) approx(1).H/diff(Z(k,:))];
      strop = 'xyz';
      for d=1:dim
        if (strcmp(phi,'rbfqr'))
          if (dim==2)
            [Bgrad{d},approx(1).Psi] = RBF_QR_diffmat_2D( ...
              strop(d),xref,approx(1).Psi);
          else
            [Bgrad{d},approx(1).Psi] = RBF_QR_diffmat_3D( ...
              strop(d),xref,approx(1).Psi);
          end  
        else
          Bgrad{d} = RBFmat(phi,ep,re,'1',d);
        end  
        %Bgrad{d} = RBFMother(phi,ep,re,'1',d, xref, xc);
        
        if (pdeg >= 0)
          deriv = zeros(1,dim); deriv(d) = 1;
          Pgrad = polyMat(xref,pdeg,deriv);
          Bgrad{d} = [Bgrad{d} Pgrad];
        end  
        Egrad{d} = (Bgrad{d}/approx(1).U)/approx(1).L; % df/dx', df/dy', df/dz' 
        Egrad{d}(:,approx(1).piv) = Egrad{d};
      end
      
      %
      % Transformation of gradients from reference to global coords
      %
%         Qtloc = T(k).Q';  
      Qtloc = T(k).Q';
      for d=1:dim
        E1{d} = zeros(length(ind),nc);
        for s = 1:dim
          E1{d} = E1{d} + dscale(s)*Qtloc(s,d)*Egrad{s}; % df/dx, df/dy, df/dz
        end
      end     
      
      %
      % Compute Hessian in reference coordinates (normalized patch)
      %
      for d_1=1:dim
      for d_2=1:dim
        if (strcmp(phi,'rbfqr'))
          if (dim==2)
            [B2{d_1,d_2},approx(1).Psi] = RBF_QR_diffmat_2D( ...
              [strop(d_1) strop(d_2)],xref,approx(1).Psi);
          else
            [B2{d_1,d_2},approx(1).Psi] = RBF_QR_diffmat_3D( ...
              [strop(d_1) strop(d_2)],xref,approx(1).Psi);
          end  
        else
          B2{d_1,d_2} = RBFmat(phi,ep,re,'1',d_1);
        end 
        %Bgrad{d} = RBFMother(phi,ep,re,'1',d, xref, xc);
        
        if (pdeg >= 0)
          deriv = zeros(dim,dim); deriv(d_1,d_2) = 1;
          P2 = polyMat(xref,pdeg,deriv);
          B2{d_1,d_2} = [B2{d_1,d_2} P2];
        end  
        E2{d_1,d_2} = (B2{d_1,d_2}/approx(1).U)/approx(1).L; % Hessian, d2f/dx'x', d2f/dx'y' etc. 
        E2{d_1,d_2}(:,approx(1).piv) = E2{d_1,d_2};
      end
      end
      
      %
      % Transformation of Hessian from reference to global coords
      %
      for d_1 = 1:dim
      for d_2 = 1:dim
        E2_global{d_1,d_2} = zeros(length(ind),nc);
        for s_1 = 1:dim
        for s_2 = 1:dim
          E2_global{d_1,d_2} = E2_global{d_1,d_2} + dscale(s_1)*dscale(s_2)*Qtloc(s_1,d_1)*Qtloc(s_2,d_2) ...
              *E2{s_1,s_2}; % Global Hessian      
        end
        end
      end
      end
      
      if (op=='L')
      %
      % Compute the Laplacian in global coords
      %
      EL = E2_global{1,1};
          for d=2:dim
            EL = EL + E2_global{d,d};
          end
      end
    end
    %
    % Form the local matrices
    %
    if (op=='0') % interpolation
        
      Eloc = spdiags(pu(k).w,0,nn,nn)*E;
      
    elseif (op=='J') % Jacobian
        Eloc{d} = zeros(length(ind),nc);
        for d = 1:dim
            Eloc{d} = spdiags(pu(k).wgrad{d},0,nn,nn)*E + spdiags(pu(k).w,0,nn,nn)*E1{d};
        end 
        
    elseif (op=='H') % Hessian
        Eloc{d_1,d_2} = zeros(length(ind),nc);
        for d_1 = 1:dim
        for d_2 = 1:dim
            Eloc{d_1,d_2} = spdiags(pu(k).wH{d_1,d_2},0,nn,nn)*E + ...
                spdiags(pu(k).wgrad{d_2},0,nn,nn)*E1{d_1} + spdiags(pu(k).wgrad{d_1},0,nn,nn)*E1{d_2} + ...
                spdiags(pu(k).w,0,nn,nn)*E2_global{d_1,d_2};
        end
        end
        
    elseif (op=='L') % Laplacian
      %
      % We are doing the Laplacian wLap*E + 2 gradw . gradE + w*ELap 
      %
      Eloc = spdiags(pu(k).w,0,nn,nn)*EL + spdiags(pu(k).wLap,0,nn,nn)*E;
      for d=1:dim
        Eloc = Eloc + 2*spdiags(pu(k).wgrad{d},0,nn,nn)*E1{d};
      end  
    end
    %
    % Add to the global matrix
    %
    ic = (nc*(k-1)+1:nc*k);
    if (op=='H')
        for d_1 = 1:dim
            for d_2 = 1:dim
                Eout{d_1,d_2}(ind,ic) = Eout{d_1,d_2}(ind,ic) + Eloc{d_1,d_2};
            end
        end
    elseif (op=='J')
        for d = 1:dim
            Eout{d}(ind,ic) = Eout{d}(ind,ic) + Eloc{d};
        end
    else
         Eout(ind,ic) = Eout(ind,ic) + Eloc;
    end
    
  end
  