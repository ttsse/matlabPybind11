%
% Point generation routine
%
function [data] = getPts(geom,N,n,C,R,mode,extCoeff)
    dim = size(C,2);
    if strcmp(geom,'ball')
        dimRat = [1 1.2*4/pi 1.2*8/(4*pi/3)];       % Ratios between rectangle/circle, cube/sphere area and volume
        dimACoeff = [1 pi (4/3)*pi];                % constant for size of domain
        dimPCoeff = [1 2*pi 4*pi];                  % constant for computing size of boundary
        h = R/(0.5*N^(1/dim) - extCoeff*n^(1/dim)); % approximate fill distance
        xB = zeros(0,dim);
        %
        % If using the unfitted method, no points on the boundary
        %
        if ~strcmp(mode,'unfitted')
            % Boundary point generation
            if dim == 1
                xB = [-R, R]';
            elseif dim == 2
                Nb = ceil((dimPCoeff(dim)*R^(dim-1))/(dimACoeff(dim-1)*h^(dim-1)));
                xB = [R*cos(linspace(-pi,pi-(2*pi)/(Nb+1),Nb)'), R*sin(linspace(-pi,pi-(2*pi)/(Nb+1),Nb)')];
            elseif dim == 3
                Nb = ceil((dimPCoeff(dim)*R^(dim-1))/(dimACoeff(dim-1)*(0.5*h)^(dim-1)));
                % Fibonacci points on sphere
                ratio = 1+sqrt(5);
                ind = [0:Nb-1]' + 0.5;
                theta = pi*ratio*ind;
                fi = acos(1-(2*ind)/(Nb));
                xB = [R.*sin(fi).*cos(theta), R.*sin(fi).*sin(theta), R.*cos(fi)];
            end
        end
        Nb = length(xB);
        % Interior point generation
        Ncube = ceil(dimRat(dim)*(N-Nb));
        x = 2*(R+n^(1/dim)*h*extCoeff)*(halton(Ncube,dim)-0.5);
        r2 = sqrt(sum(x.^2,2));
        pos = find(r2<=(R+n^(1/dim)*h*extCoeff));
        pos = pos(1:N-Nb);
        x = x(pos,:) + C;
        %
        % Organize outputs including labels for points inside, outside and on the boundary
        %
        data.nodes = [x; xB];
        data.inner = find(sqrt(sum((x+C).^2,2))<=R);
        data.outer = setdiff(1:size(x,1),data.inner)';
        data.bnd = [size(x,1) + 1:size(x,1) + size(xB,1)]';
        %
        % Domain volume and area measures used for the scaling of the LS problem
        %
        data.Vol = (dimACoeff(dim)*R^dim);
        data.Area = (dimPCoeff(dim)*R^(dim-1));
    elseif strcmp(geom,'cube')
        dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];          % constants for getting cube side length /2
        h = R/(0.5*N^(1/dim) - extCoeff*n^(1/dim));   % approximate fill distance
        xB = zeros(0,dim);
        %
        % If using the unfitted method, no points on the boundary
        %
        if ~strcmp(mode,'unfitted')
            % Boundary point generation
            if dim == 1
                xB = [-R, R]';
            elseif dim == 2
                Nb = ceil(((2*dimLCoeff(dim)*R)^(dim-1))/(h^(dim-1)));
                XYZLim = [-dimLCoeff(dim)*R dimLCoeff(dim)*R];
                xB = [linspace(XYZLim(1),XYZLim(2),Nb)', ones(Nb,1)*XYZLim(1);...
                      ones(Nb,1)*XYZLim(1), linspace(XYZLim(1),XYZLim(2),Nb)';...
                      linspace(XYZLim(1),XYZLim(2),Nb)', ones(Nb,1)*XYZLim(2);...
                      ones(Nb,1)*XYZLim(2), linspace(XYZLim(1),XYZLim(2),Nb)'];
                xB = unique(xB,'rows');
                xB = xB + C;
            elseif dim == 3
                NbF = ceil(((2*dimLCoeff(dim)*R)^(dim-1))/(h^(dim-1)));     % Number of points on each face
                NbE = ceil(2*dimLCoeff(dim)/h);                             % Number of points on each edge
                XYZLim = [-dimLCoeff(dim)*R dimLCoeff(dim)*R];
                xB = [linspace(XYZLim(1),XYZLim(2),NbE)', ones(NbE,1)*XYZLim(1), ones(NbE,1)*XYZLim(1);...
                      ones(NbE,1)*XYZLim(1), linspace(XYZLim(1),XYZLim(2),NbE)', ones(NbE,1)*XYZLim(1);...
                      linspace(XYZLim(1),XYZLim(2),NbE)', ones(NbE,1)*XYZLim(2), ones(NbE,1)*XYZLim(1);...
                      ones(NbE,1)*XYZLim(2), linspace(XYZLim(1),XYZLim(2),NbE)', ones(NbE,1)*XYZLim(1)];
                xB = [xB; [1,1,-1].*xB];
                Ry = eye(dim); Ry(1,1) = cos(pi/2); Ry(dim,dim) = cos(pi/2); Ry(dim,1) = -sin(pi/2); Ry(1,dim) = sin(pi/2); % Rotate 90 degrees along y-axis
                xB = [xB; (Ry*xB')'];
                xB = uniquetol(xB,1e-14,'ByRows',true);
                NbE = size(xB,1);
                Rx = eye(dim); Rx(dim-1,dim-1) = cos(pi/2); Rx(dim,dim) = cos(pi/2); Rx(dim-1,dim) = -sin(pi/2); Rx(dim,dim-1) = sin(pi/2); % Rotate 90 degrees
                xBF = [2*(dimLCoeff(dim)*R)*(halton(NbF,dim-1)-0.5), ones(NbF,1)*XYZLim(1)];
                xBF = [xBF; [1,1,-1].*xBF];
                xB = [xB; xBF];
                xB = [xB; (Rx*xBF')'];
                xB = [xB; (Ry*xBF')'];
            end
        end
        Nb = length(xB);
        % Interior point generation
        x = 2*(dimLCoeff(dim)*R+dimLCoeff(dim)*n^(1/dim)*h*extCoeff)*(halton(N-Nb,dim)-0.5);
        r2 = sqrt(sum(x.^2,2));
        x = x + C;
        %
        % Organize outputs including labels for points inside, outside and on the boundary
        %
        data.nodes = [x; xB];
        data.inner = find(all(abs(x-C)<=dimLCoeff(dim)*R,2));
        data.outer = setdiff(1:size(x,1),data.inner)';
        data.bnd = [size(x,1) + 1:size(x,1) + size(xB,1)]';
        %
        % Domain volume and area measures used for the scaling of the LS problem
        %
        data.Vol = (2*dimLCoeff(dim)*R)^dim;
        data.Area = 2*dim*(2*dimLCoeff(dim)*R)^(dim-1);
    end
end