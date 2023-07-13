% Evaluate Wendland's C� function and its derivatives
% Modified for cylindrical patches 2019
% The bump function added 2019

%  (C) Alfa Heryudono 2010, Victor Shcherbakov, Elisabeth Larsson 2015, 2018

function [phi,phigrad,phiLap,phiH] = cylWeightRBF(psi,C,R,Z,T,xloc,plist)
%
% First compute all derivatives on a unit domain. Then scale them down
% according to Z and R. We still need to make locpts to be unit.
%

switch psi
    case 'bump'
        funphi  = @(r) exp(-1./(1-r.^2)).*(r<1);
        % This is the gradient divided by r
        gradphi = @(r) 2./(1-r.^2).^2.*funphi(r);
        phirr   = @(r) (2*(1-r.^2)+8*r.^2)./(1-r.^2).^3.*funphi(r);
    case 'wendland_c2'
        funphi = @(r) (1-r).^4.*(4*r + 1) .* (r<1);
        gradphi = @(r) (4*(r - 1).^4 + 4*(4*r + 1).*(r - 1).^3 ).* (r<1);
        phirr = @(r) (32*(r - 1).^3 + 12*(4*r + 1).*(r - 1).^2).* (r<1);
    otherwise
        error(['Weight function ' psi ' is not implented'])
end

  [npt,np] = size(plist);
  dim = size(C,2);

  empty = spalloc(npt,np,nnz(plist));  
  [phi,phiLap] = deal(empty);
  for d=1:dim
    phigrad{d} = empty;
    for k=1:dim
      phiH{k,d} = empty;
    end  
  end
  %
  % Go through each patch
  %
  for i = 1:np
    ind = find(plist(:,i)); % Points in this patch
    num = sum(plist(ind,:),2); % Patches this point are in
    only = ind(find(num==1)); % Weight one in this patch
    lind = length(ind);
    szx = size(xloc(i).pts);
    pos = find(num > 1);
    xloc(i).pts = xloc(i).pts(find(num > 1),:); 
    ind  = ind(find(num > 1)); % Local overlap points
    
    phi(only,i) = 1; % All derivatives are zero in this case
    %
    % If there are points in this patch to compute for
    %
    H = 2; % Unit cylinder height used for defining weight function
    if (length(ind)>0)
      locgrad = zeros(length(ind),dim);
      %
      % The local coordinate system is around the horizontal patch center.
      % The vertical direction has an offset. We compute in the actual patch
      % while the functions are defined on r<1 and |z|<1.
      %
      dr = xloc(i).pts(:,1:dim-1)/R(i); % One or two columns
      dz = (xloc(i).pts(:,dim)-mean(Z(i,:)))/diff(Z(i,:))*H;
      
      r = sqrt(sum(dr.^2,2));
      pos = find(r>=1);
      r(pos) = 2; % To ensure we get zero. Any number witn an eps
                  % distance from one works
      
      rz = abs(dz);
      pos = find(rz>=1);
      rz(pos) = 2; % To ensure we get zero. Any number witn an eps
                  % distance from one works
      
      phi(ind,i) = funphi(r).*funphi(rz);

      locgrad(:,1:dim-1) = dr.*gradphi(r).*funphi(rz)/R(i);
      locgrad(:,dim) = dz.*gradphi(rz).*funphi(r)/diff(Z(i,:))*H;
          
      grad = locgrad*T(i).Q';
      for d=1:dim
        phigrad{d}(ind,i) = grad(:,d);
      end  
      %
      % The Laplacian is Delta_xy*phi(z)+Delta_z*phi(r)
      % The radial part becomes 1/r d/dr(r dphi/dr) = 1/r phir + phirr
      % 
      phiLap(ind,i) = ((dim-2)*gradphi(r) + phirr(r)).*funphi(rz)/R(i)^2;
      %if (dim==2)
      %  phiLap(ind,i) = phirr(r).*funphi(rz)/R(i)^2;
        %phiLap(ind,i) = phiLap(ind,i) + ...
        %    phirr(rz).*funphi(r)/diff(Z(i,:)).^2*H^2;
        %else
        % Radial Laplacian in two coordinates
        %phiLap(ind,i) = (gradphi(r) + phirr(r)).*funphi(rz)/R(i)^2;
        %phiLap(ind,i) = phiLap(ind,i) + ...
        %    phirr(rz).*funphi(r)/diff(Z(i,:)).^2*H^2; 
        %end
      phiLap(ind,i) = phiLap(ind,i) + ...
            phirr(rz).*funphi(r)/diff(Z(i,:)).^2*H^2; 
      %
      % The Hessian elements in the radial directions are 
      % (1-x_i^2/r^2)*phigrad/r + phirr*x_i^2/r^2 and
      % -x_ix_j/r^2*phigrad/r+phirr*x_i*y_i/r^2
      % Scaling to be included
      %
      % Mixed radial and z-directions
      % x_i phigrad/r* z *phigrad/r Can be gotten from the grads.
      %
      % Second derivative in z
      % phirr(rz).*funphi(r)/diff(Z(i,:)).^2*H^2
      coeff = eye(dim-1);
      ratio = dr./r;
      locH = zeros(size(r,1),dim,dim);
      for d=1:dim-1
        for k=1:dim-1
          locH(:,d,k) = ((coeff(d,k)-ratio(:,d).*ratio(:,k)) .* ...
              gradphi(r)+phirr(r).*ratio(:,d).*ratio(:,k)).* ...
              funphi(rz)/R(i).^2;   
        end
        locH(:,d,dim) = dr(:,d).*gradphi(r).* ...
            dz.*gradphi(rz)/R(i)/diff(Z(i,:))*H;
        locH(:,dim,d) = locH(:,d,dim);
      end
      locH(:,dim,dim) = phirr(rz).*funphi(r)/diff(Z(i,:)).^2*H^2; 
      %
      % Now we need to apply QHQ' to get the Hessian in physical space
      % We need to think carefully so we multiply in the right way
      %
      if (size(locH,1)==1)
        tmp = squeeze(locH(1,:,:));
        locH = T(i).Q*tmp*T(i).Q';
        locH(1,1:dim,1:dim) = tmp;
        locH = locH(1,1:dim,1:dim);
      else % More elements  
        for d=1:dim % Multiply each row with Q' from the right
          locH(:,d,:) = squeeze(locH(:,d,:))*T(i).Q';
        end
        for d=1:dim % Multiply each column with Q from the left
          %locH(:,:,d) = (T(i).Q*locH(:,:,d)')';
          locH(:,:,d) = locH(:,:,d)*T(i).Q';
        end
      end  
      %
      % Assign the output variable
      %
      for d=1:dim
        for k=1:dim
            phiH{k,d}(ind,i) = locH(:,k,d);
        end
      end  
    end
    %
    % Have validated that trH == Laplacian before and after transformation
    %   
  end










