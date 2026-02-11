% -------------------------------------------------------------------------
% plotPtch.m -- Plotting patches in 1D, 2D or 3D.
% Inputs     -- ptch    -> Patch structure. See getPtch.
%               geom    -> Geometry to solve, geom = "cube" or "ball".
%               C       -> Double 1xd array of domain centre in 
%                                      d dimensions.
%               R       -> Double. Domain radius.
% Syntax     -- pars = setPars;
%               n = nchoosek(pars.rbfdeg+pars.dim,pars.dim);
%               [ptch, dataX, dataY] = getPtch([0,0],1,n,pars);
%               plotPtch(ptch,'ball',[0,0],1);
%
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function [] = plotPtch(ptch,geom,C,R)
    theta = linspace(0,2*pi,100)';
    dim = size(C,2);
    figure()
    %
    % Plot patches
    %
    if dim == 1
        altern = 0;
        for i = 1:size(ptch.C,1)
            plot(ptch.C(i,1),altern,'ro');
            hold on
            x = [ptch.C(i,1) + ptch.R(i); ptch.C(i,:) - ptch.R(i)];
            plot(x(:,1),altern.*ones(2,1),'k-')
            altern = mod(i,2)*1e-2;
        end
    elseif dim == 2
        for i = 1:size(ptch.C,1)
            plot(ptch.C(i,1),ptch.C(i,2),'ro');
            hold on
            x = ptch.C(i,:) + [ptch.R(i)*cos(theta), ptch.R(i)*sin(theta)];
            plot(x(:,1),x(:,2),'k-')
        end
    elseif dim == 3
        [x,y,z] = sphere(10);
        for i = 1:size(ptch.C,1)
            plot3(ptch.C(i,1),ptch.C(i,2),ptch.C(i,3),'ro');
            hold on
            surf(x.*ptch.R(i)+ptch.C(i,1),y.*ptch.R(i)+ptch.C(i,2),z.*ptch.R(i)+ptch.C(i,3),...
                'FaceColor','none');
        end
    end
    %
    % Plot geometry frame
    %
    if strcmp(geom,"cube")
        dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];   
        a = dimLCoeff(dim)*R;
        if dim == 1
            x = [C + a; ...
                 C - a];
            plot(x(:,1),zeros(size(x,1)),'ko');
        elseif dim==2
            x = [C + [a a]; ...
                 C + [a -a]; ...
                 C + [-a -a]; ...
                 C + [-a a]; ...
                 C + [a a]];
            plot(x(:,1),x(:,2),'k-');
        elseif dim==3
            x = [C + [a a a]; ...
                 C + [a -a a]; ...
                 C + [-a -a a]; ...
                 C + [-a a a]; ...
                 C + [a a a]];
            x = [x; [x(:,1), x(:,2), x(:,3)-2*a]];
            plot3(x(:,1),x(:,2),x(:,3),'k-','LineWidth',2);
            for i = 1:size(x,1)/2
                xPlot = [x(i,:); x(i+size(x,1)/2,:)];
                plot3(xPlot(:,1),xPlot(:,2),xPlot(:,3),'k-','LineWidth',2);
            end
        end
    elseif strcmp(geom,"ball")
        if dim == 1
            x = [C + R; ...
                 C - R];
            plot(x(:,1),zeros(size(x,1)),'ko');
        elseif dim == 2
            theta = linspace(0,2*pi,1000);
            plot(R.*cos(theta),R.*sin(theta),'k-')
        elseif dim == 3
            [x,y,z] = sphere(20);
            surf(x.*R+C(1),y.*R+C(2),z.*R+C(3),'FaceColor','none','LineWidth',2);
        end
    end
    axis equal
end