% -------------------------------------------------------------------------
% computeDifferror.m -- Test function which computes the L2 and Linf error
%                       of numerical interpolation operators and derivatives
%                       Test function used f = sin(pi*x*y*z).
% Copyright (c) 2023 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------

function differror = computeDifferror(xc,xi,xb,Ei,Eb,Jb,Hi)

%         fscale = [194.6542 136.2784 925.8762];
        fscale = [abs(max(xb(:,1))-min(xb(:,1))) ...
            abs(max(xb(:,2))-min(xb(:,2))) abs(max(xb(:,3))-min(xb(:,3)))];
        syms x1 y1 z1

        u_test = sin(pi*((x1-mean(xb(:,1)))./fscale(1))*((y1-mean(xb(:,2)))./fscale(2))*((z1-mean(xb(:,3)))./fscale(3)));
        % u_test = (x1 + y1 + z1).^3;
        
        fun = matlabFunction(u_test);

        dx_fun = matlabFunction(diff(u_test,x1));
        dy_fun = matlabFunction(diff(u_test,y1));
        dz_fun = matlabFunction(diff(u_test,z1));

        dxx_fun = matlabFunction(diff(u_test,x1,x1));
        dyy_fun = matlabFunction(diff(u_test,y1,y1));
        dzz_fun = matlabFunction(diff(u_test,z1,z1));

        dxy_fun = matlabFunction(diff(u_test,x1,y1));
        dyz_fun = matlabFunction(diff(u_test,y1,z1));
        dxz_fun = matlabFunction(diff(u_test,x1,z1));

        dyx_fun = matlabFunction(diff(u_test,y1,x1));
        dzy_fun = matlabFunction(diff(u_test,z1,y1));
        dzx_fun = matlabFunction(diff(u_test,z1,x1));
                
        differror.i = norm(Ei*fun(xc(:,1), xc(:,2), xc(:,3)) - fun(xi(:,1),xi(:,2),xi(:,3)),2)/norm(fun(xi(:,1),xi(:,2),xi(:,3)));
        differror.b = norm(Eb*fun(xc(:,1), xc(:,2), xc(:,3)) - fun(xb(:,1),xb(:,2),xb(:,3)),2)/norm(fun(xb(:,1),xb(:,2),xb(:,3)));

        differror.dxb = norm(Jb{1}*fun(xc(:,1), xc(:,2), xc(:,3)) - dx_fun(xb(:,1),xb(:,2),xb(:,3)),2)/norm(dx_fun(xb(:,1),xb(:,2),xb(:,3)));
        differror.dyb = norm(Jb{2}*fun(xc(:,1), xc(:,2), xc(:,3)) - dy_fun(xb(:,1),xb(:,2),xb(:,3)),2)/norm(dy_fun(xb(:,1),xb(:,2),xb(:,3)));
        differror.dzb = norm(Jb{3}*fun(xc(:,1), xc(:,2), xc(:,3)) - dz_fun(xb(:,1),xb(:,2),xb(:,3)),2)/norm(dz_fun(xb(:,1),xb(:,2),xb(:,3)));

        differror.dxx = norm(Hi{1,1}*fun(xc(:,1), xc(:,2), xc(:,3)) - dxx_fun(xi(:,1),xi(:,2),xi(:,3)),2)/norm(dxx_fun(xi(:,1),xi(:,2),xi(:,3)));
        differror.dyy = norm(Hi{2,2}*fun(xc(:,1), xc(:,2), xc(:,3)) - dyy_fun(xi(:,1),xi(:,2),xi(:,3)),2)/norm(dyy_fun(xi(:,1),xi(:,2),xi(:,3)));
        differror.dzz = norm(Hi{3,3}*fun(xc(:,1), xc(:,2), xc(:,3)) - dzz_fun(xi(:,1),xi(:,2),xi(:,3)),2)/norm(dzz_fun(xi(:,1),xi(:,2),xi(:,3)));

        differror.dxy = norm(Hi{1,2}*fun(xc(:,1), xc(:,2), xc(:,3)) - dxy_fun(xi(:,1),xi(:,2),xi(:,3)),2)/norm(dxy_fun(xi(:,1),xi(:,2),xi(:,3)));
        differror.dyz = norm(Hi{2,3}*fun(xc(:,1), xc(:,2), xc(:,3)) - dyz_fun(xi(:,1),xi(:,2),xi(:,3)),2)/norm(dyz_fun(xi(:,1),xi(:,2),xi(:,3)));
        differror.dxz = norm(Hi{1,3}*fun(xc(:,1), xc(:,2), xc(:,3)) - dxz_fun(xi(:,1),xi(:,2),xi(:,3)),2)/norm(dxz_fun(xi(:,1),xi(:,2),xi(:,3)));

        differror.dyx = norm(Hi{2,1}*fun(xc(:,1), xc(:,2), xc(:,3)) - dyx_fun(xi(:,1),xi(:,2),xi(:,3)),2)/norm(dyx_fun(xi(:,1),xi(:,2),xi(:,3)));
        differror.dzy = norm(Hi{3,2}*fun(xc(:,1), xc(:,2), xc(:,3)) - dzy_fun(xi(:,1),xi(:,2),xi(:,3)),2)/norm(dzy_fun(xi(:,1),xi(:,2),xi(:,3)));
        differror.dzx = norm(Hi{3,1}*fun(xc(:,1), xc(:,2), xc(:,3)) - dzx_fun(xi(:,1),xi(:,2),xi(:,3)),2)/norm(dzx_fun(xi(:,1),xi(:,2),xi(:,3)));
end