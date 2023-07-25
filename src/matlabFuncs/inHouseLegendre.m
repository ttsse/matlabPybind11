% ----------------------------------------------------------------------------------
% inHouseLegendre.m -- Computing associated Legendre polynomials of degree l and
%                      order m = 0,..,length(x). The first and second derivative are
%                      also computed (dPdth, dPdth2), as well as combinations A, B.
%                      Option to compute normalised values exists.
% Copyright (c) 2023 Andreas Michael <andreas.michael@it.uu.se>, 
%               Elisabeth Larsson <elisabeth.larsson@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% ----------------------------------------------------------------------------------

function [P dPdth pOverS dPdth2 A B] = inHouseLegendre(l, x, normal)

% Refer to [1], [2] and [3] for details on formulas used to compute the polynomials.
% 
% References:
%     [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%         Functions", Dover Publications, 1965, Ch. 8.
%     [2] "Associated Legendre Polynomials." Wikipedia, Wikimedia Foundation, 
%         23 Nov. 2022, en.wikipedia.org/wiki/Associated_Legendre_polynomials. 
%         Accessed 25 Jul. 2023.
%     [3] D. B. Westra, "Identities and properties for associated Legendre functions",
%         www.mat.univie.ac.at/~westra/associatedlegendrefunctions.pdf. 
%         Accessed 25 Jul 2023.

if max(x) > 1 || min(x) < -1
    error("Input vector x with values between -1 and 1!")
end

if l == 0
    error("Only works with integer values of l > 0")
end

x = reshape(x,[1, length(x)]);

P = zeros(l+1, length(x));
pOverS = zeros(l+1, length(x));
pOverS2 = zeros(l+1, length(x));
dPdth = zeros(l+1, length(x));
dPdth2 = zeros(l+1, length(x));
A = zeros(l+1, length(x));
B = zeros(l+1, length(x));

P(l+1, :) = (factorial(2*l)/((2^l)*factorial(l)))*((-1)^l)*((1-x.^2).^(l/2));

idx = find(abs(x)~=1);
ids = setdiff(1:length(x),idx);

P(l, idx) = -(x(idx)./((1-x(idx).^2).^(0.5))).*P(l+1, idx);

for m = l-2:-1:0
    C = (1/((l+m+1)*(l-m)));
    P(m+1, idx) = -C.*((2*(m+1).*x(idx).*(1-x(idx).^2).^(-0.5)).*P(m+2,idx) + P(m+3,idx));
end

P(1, ids) = x(ids).^l;

pOverS(:, idx) = P(:,idx)./sqrt(1 - x(idx).^2);
pOverS(2, ids) = -(x(ids).^(l-1))*l*(l+1)*.5;
pOverS2(:, idx) = P(:,idx)./(1 - x(idx).^2);
pOverS2(3, ids) = (x(ids).^(l-2))*(l-1)*l*(l+1)*(l+2)*.125;
pOverS2(2, idx) = 0;
pOverS2(1, idx) = 0;

if nargin > 2 && normal == "unnorm"

    dPdth(1,:) = -P(2,:);
    dPdth(l+1,:) = l.*P(l,:);
    
    for m = 1:l-1
        C = (l+m)*(l-m+1);
        dPdth(m+1,:) = 0.5.*(C.*P(m,:) - P(m+2,:));
    end
    
    dPdth2(1,:) = -dPdth(2,:);
    dPdth2(l+1,:) = l.*dPdth(l,:);
    
    for m = 1:l-1
        C = (l+m)*(l-m+1);
        dPdth2(m+1,:) = 0.5.*(C.*dPdth(m,:) - dPdth(m+2,:));
    end

elseif nargin > 2 && normal == "norm"
    
    for m = l:-1:0
        coeff = ((-1).^m)*sqrt(((2*l + 1)/2)*(factorial(l-m)/factorial(l+m)));
        P(m+1,:) = coeff*P(m+1,:); 
        pOverS(m+1,:) = coeff*pOverS(m+1,:);
        pOverS2(m+1,:) = coeff*pOverS2(m+1,:);
    end

    dPdth(1,:) = sqrt(l*(l+1)).*P(2,:);
    dPdth(l+1,:) = -sqrt(l/2).*P(l,:);
    for m = 1:l-1
        C1 = sqrt((l+m)*(l-m+1));
        C2 = sqrt((l-m)*(l+m+1));
        dPdth(m+1,:) = -0.5.*(C1.*P(m,:) - C2.*P(m+2,:));
    end

    dPdth2(1,:) = sqrt(l*(l+1)).*dPdth(2,:);
    dPdth2(l+1,:) = -sqrt(l/2).*dPdth(l,:);
    for m = 1:l-1
        C1 = sqrt((l+m)*(l-m+1));
        C2 = sqrt((l-m)*(l+m+1));
        dPdth2(m+1,:) = -0.5.*(C1.*dPdth(m,:) - C2.*dPdth(m+2,:));
    end
   
    if (l==1)
      A(l+1,:) = -P(2,:);
      B(l+1,:) = -P(2,:);
    else
      A(l+1,:) = -(l^2)*pOverS2(l+1,:) ...
          +x*sqrt(0.5*l).*pOverS(l,:);
      B(l+1,:) = -l*(pOverS2(l+1,:) ...
                      -x.*sqrt(0.5*l).*pOverS(l,:));
    end
    A(1,:) = -x.*sqrt(l*(l+1)).*pOverS(2,:);

    B(1,:) = 0;
    for m=1:l-1
      A(m+1,:) = -m*(m-1)*pOverS2(m+1,:)-m*P(m+1,:) ...
          -sqrt(l+m+1)*sqrt(l-m)*x.*pOverS(m+2,:);
      
      B(m+1,:) = m*(m-1)*pOverS2(m+1,:)-m*m*P(m+1,:) ...
          -m*sqrt(l+m+1)*sqrt(l-m)*x.*pOverS(m+2,:);
    end 
else
    error("Normalisation argument not recognised")
end

end