% -------------------------------------------------------------------------
% RBFsolver.m -- Numerically compute the solution to a Poisson problem using 
%                RBF based methods. Modules cpp, np, sp required (See
%                importModules) if pars.cppOn = 1.
% Inputs         -- pars    -> Parameter structure, see setPars.
%                   cpp     -> Shared python library generated using
%                              pybind11. See importModules.
%                   np      -> Numpy library. See importModules.
%                   sp      -> Scipy library. See importModules.
% Outputs        -- l2Error -> Double. Discrete numerical error.
%                   h       -> Double. Fill distance between node points.
% Syntax         -- pathSet;
%                   pars = setPars;
%                   [cpp, np, sp] = importModules;
%                   [l2Error, h] = RBFsolver(pars,varargin)   
%
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function [l2Error, h] = RBFsolver(pars,varargin)
%
% Ensure relevant libraries are provided
%
if pars.cppOn == 1 
    if nargin == 4 && (strcmp(underlyingType(varargin{1}),'py.module') && strcmp(underlyingType(varargin{2}),'py.module') && strcmp(underlyingType(varargin{3}),'py.module')) 
        [cpp, np, sp] = varargin{1:3};
    else
        error("RBFsolver: If parameter cppOn = 1, function requires python modules cpp, np, sp. Run [cpp, np, sp] = importModules and then RBFsolver(pars,cpp,np,sp).");
    end
end
% 
% Point number for RBF-FD method choice
%
if strcmp(pars.method,'FD')
    if pars.pdeg == -1
        n = nchoosek(pars.rbfdeg+pars.dim-1,pars.dim); % Stencil size in smooth RBF case
    else
        n = 2*nchoosek(pars.pdeg+pars.dim,pars.dim);   % Stencil size in phs + poly case
    end
    if pars.N < n && strcmp(pars.method,'FD')
        warning("RBFFD: The number of centre points are less than the stencil points required for the given polynomial order. N is increased.")
        pars.N = n;
    end
    extCoeff = pars.extCoeff*(strcmp(pars.mode,"unfitted"));
elseif strcmp(pars.method,'PUM')
    n = nchoosek(pars.rbfdeg+pars.dim,pars.dim);
end
%
% Place N centre points and M evaluation points in geom with centre C and radius R
%
C = zeros(1,pars.dim);
if strcmp(pars.geom,'cube')
    R = 1*(pars.dim)^(1/2);
elseif strcmp(pars.geom,'ball')
    R = 1;
else
    error("RBFsolver: Requested geometry not implemented yet");
end
if strcmp(pars.method,'FD')
    %
    % Get centre points (xc)
    %
    [dataX] = getPts(pars.geom,pars.N,n,C,R,pars.mode,extCoeff);
    xc = dataX.nodes;                                     
    % 
    % Get evaluation points (xe)
    %
    if strcmp(pars.mode,"collocation")       
        M = pars.N;                             
        dataY = dataX;
        xe = dataY.nodes;
    else
        M = ceil(pars.q*pars.N);                              
        [dataY] = getPts(pars.geom,M,n,C,R,"fitted",0);
        %
        % Move evaluation points inside (and on the boundary) to the closest center point.
        % Make sure not to move two evaluation points to the same center point. All
        % center points should be moved to.
        %
        if pars.mvCentres && ~strcmp(pars.mode,"collocation")
            dataY = movePts(dataX,dataY);
        end
        xe = dataY.nodes;
    end
    %
    % Ensure center points are not too far from boundary and make evaluation point - stencil list 
    %
    if strcmp(pars.mode,"unfitted")
        xc = xc(unique(knnsearch(xc,xe,'K',ceil(extCoeff*n+pars.ep))),:);
        N = size(xc,1);
    else 
        N = pars.N;
    end
    ptStencilList = knnsearch(xc,xe,'K',1);
elseif strcmp(pars.method,'PUM')
    %
    % Get patches 
    %
    ptch = getPtch(pars.geom,pars.P,C,R,pars.del);
    P = length(ptch.R);
    %
    % Get centre points (xc)
    %
    if strcmp(pars.mode,"unfitted")
        dataX.nodes = [];
        for i = 1:P
            [ptch.xc(i)] = getPts("ball",n,0,ptch.C(i,:),ptch.R(i),"unfitted",0);
            xCglobalId{i} = [(i-1)*n+1:i*n]';
            dataX.nodes = [dataX.nodes; ptch.xc(i).nodes];
        end
        xc = dataX.nodes;
        [ptch.xc.globalId] = xCglobalId{:};
        if strcmp(pars.geom,"ball")
            dataX.inner = find(sqrt(sum((xc-C).^2,2))<=R);
        elseif strcmp(pars.geom,"cube")
            dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];  
            dataX.inner = find(all(abs(xc-C)<=dimLCoeff(pars.dim)*R,2));
        end
        dataX.outer = setdiff(1:size(xc,1),dataX.inner);
        dataX.bnd = [];
    elseif strcmp(pars.mode,"collocation")
        dataX = getPts(pars.geom,n*P,0,C,R,"fitted",0);
        xc = dataX.nodes;
        for i = 1:P
            ptch.xc(i).globalId = find(sqrt(sum((xc - ptch.C(i,:)).^2,2)) <= ptch.R(i));
            ptch.xc(i).nodes = xc(ptch.xc(i).globalId,:);
        end
    elseif strcmp(pars.mode,"fitted") 
        dataX = getPts(pars.geom,n*P,0,C,R,"fitted",0);
        xc = dataX.nodes;
        for i = 1:P
            ptch.xc(i).globalId = find(sqrt(sum((xc - ptch.C(i,:)).^2,2)) <= ptch.R(i));
            ptch.xc(i).nodes = xc(ptch.xc(i).globalId,:);
        end
    end
    %
    % Get evaluation points (xe)
    %
    q = max(pars.q*double(~strcmp(pars.mode,"collocation")),1);
    M = ceil(n*P*q);  
    if ~strcmp(pars.mode,"collocation")
        dataY = getPts(pars.geom,M,0,C,R,"fitted",0);
        %
        % Move evaluation points inside (and on the boundary) to the closest center point
        %
        if pars.mvCentres && ~strcmp(pars.mode,"collocation")
            dataY = movePts(dataX,dataY);
        end
        xe = dataY.nodes;
        for i = 1:P
            ptch.xe(i).globalId = find(sqrt(sum((xe - ptch.C(i,:)).^2,2)) <= ptch.R(i));
            ptch.xe(i).nodes = xe(ptch.xe(i).globalId,:);
        end
    else
        dataY = dataX;
        xe = xc; 
        ptch.xe = ptch.xc;
    end
    N = size(xc,1);
end
%
% Constructing global approximation to evaluation and Laplace operators
%
if strcmp(pars.method,'FD')
    [L,B, Eglobal, Lglobal] = conGlobMat(dataY,M,N,n,pars.prob,pars.phi,pars.psi,pars.ep,pars.pdeg,ptStencilList,xc,xe);
elseif strcmp(pars.method,'PUM')
    [L,B, Eglobal, Lglobal] = conGlobMat(dataY,M,N,n,pars.prob,pars.phi,pars.psi,pars.ep,pars.pdeg,ptch);
end
%
% Manufactured solution to construct forcing and BC
%
if pars.dim == 1
    fun = @(x) sin(2.*pi.*x);
    lapFun = @(x) -4.*(pi.^2).*sin(2.*pi.*x);
    F = [lapFun(xe(dataY.inner,1)); fun(xe(dataY.bnd,1))];
    uExact = fun(xe(:,1));
    ucExact = fun(xc(:,1));
    lapAnalytic = lapFun(xe(dataY.inner,1));
    bndAnalytic = [0; 0];
elseif pars.dim == 2
    fun = @(x,y) sin(2.*pi.*x.*y);
    lapFun = @(x,y) - 4.*(x.^2).*(pi.^2).*sin(2.*pi.*x.*y) - 4.*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y);
    F = [lapFun(xe(dataY.inner,1),xe(dataY.inner,2)); fun(xe(dataY.bnd,1),xe(dataY.bnd,2))];
    uExact = fun(xe(:,1),xe(:,2));
    ucExact = fun(xc(:,1),xc(:,2));
    lapAnalytic = lapFun(xe(dataY.inner,1),xe(dataY.inner,2));
    bndAnalytic = fun(xe(dataY.bnd,1),xe(dataY.bnd,2));
elseif pars.dim == 3
    fun = @(x,y,z) sin(2.*pi.*x.*y.*z);
    lapFun = @(x,y,z) - 4.*(x.^2).*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(y.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(x.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z);
    F = [lapFun(xe(dataY.inner,1),xe(dataY.inner,2),xe(dataY.inner,3)); fun(xe(dataY.bnd,1),xe(dataY.bnd,2),xe(dataY.bnd,3))];
    uExact = fun(xe(:,1),xe(:,2),xe(:,3));
    ucExact = fun(xc(:,1),xc(:,2),xc(:,3));
    lapAnalytic = lapFun(xe(dataY.inner,1),xe(dataY.inner,2),xe(dataY.inner,3));
    bndAnalytic = fun(xe(dataY.bnd,1),xe(dataY.bnd,2),xe(dataY.bnd,3));
end 
%
% LS problem scaling
%
[~,dist] = knnsearch(xc,xc,'K',2);
h = max(dist(:,2));
if pars.scaling && (strcmp(pars.mode,"unfitted") || (strcmp(pars.mode,"fitted") && strcmp(pars.bcMode,"weak")))
    Lscale = sqrt(dataY.Vol/length(dataY.inner));
    Bscale = sqrt(dataY.Area/length(dataY.bnd))*(h.^-1.5);
    L = Lscale.*L; B = Bscale.*B;
    F(dataY.inner) = Lscale.*F(dataY.inner);
    F(dataY.bnd) = Bscale.*F(dataY.bnd);
end

if strcmp(pars.mode,"fitted") && strcmp(pars.bcMode,"strong")
    Fmod = Lglobal(:,dataX.bnd)*ucExact(dataX.bnd);
    F = F-Fmod;
    L = Lglobal(:,dataX.inner);
    B = zeros(0,N);
end
%
% Solution on centre points, evaluated on Y set
%
A = [L; B];
if pars.cppOn
    testMemory(nnz(A),M,pars.memTol);
    [BI, BJ, BV] = find(A);
    row  = np.array(BI - 1).reshape(int32(-1));
    col  = np.array(BJ - 1).reshape(int32(-1));
    data = np.array(BV,pyargs('dtype', 'float64')).reshape(int32(-1));
    pyB = sp.csc_matrix({data,{row, col}}, {int32(size(A,1)), int32(size(A,2))});
    u = double(cpp.sparseSolve(pyB,np.array(F,pyargs('dtype', 'float64')),pars.debug))';
    if pars.debug && pars.display
        disp(['RBFsolver: c++ sparse QR solver solve time was ', num2str(u(end)), 's']);
        sTimeCpp = u(end);
        u(end) = [];
    end
else
    if pars.debug
        tic
        u = A\F;
        sTimeMat = toc;
    else
        u = A\F;
    end
end
%
% Fix operators to compute error measures
%
if strcmp(pars.mode,"fitted") && strcmp(pars.bcMode,"strong")
    u = [u; ucExact(dataX.bnd)];
end
ue = Eglobal*u;
l2Error = norm(abs(ue-uExact),2)/norm(uExact,2);
%
% Displaying output 
%
if pars.display
    if strcmp(pars.method,'FD')
        figure()
        hold on;
    elseif strcmp(pars.method, 'PUM')
        plotPtch(ptch,pars.geom,C,R)
        hold on
    end
    if pars.dim == 1
        evalPtPlot = plot(xe,zeros(size(xe,1),1),'b.');
        centerPtPlot = plot(xc,zeros(size(xc,1),1),'rx');
    elseif pars.dim == 2
        centerPtPlot = plot(xc(:,1),xc(:,2),'rx','MarkerSize',8,'LineWidth',3);
        if ~strcmp(pars.mode,"collocation")
            evalPtPlot = plot(xe(:,1),xe(:,2),'b.','MarkerSize',8);
        end
    elseif pars.dim == 3
        evalPtPlot = plot3(xe(:,1),xe(:,2),xe(:,3),'b.');
        centerPtPlot = plot3(xc(:,1),xc(:,2),xc(:,3),'rx');
    end
    ax = gca;
    ax.LineWidth = 1.5;
    ax.FontSize = 18;
    xlabel("x",'Interpreter','latex','FontSize',20)
    ylabel("y",'Interpreter','latex','FontSize',20,'Rotation',0)
    axis equal
    xlim([-1.65 1.65]);
    ylim([-1.5 1.3]);
    if strcmp(pars.mode,"unfitted")
        if strcmp(pars.method,'FD')
            l = legend('Node points','Eval points','Interpreter','latex','location','south','Orientation','Horizontal');
        elseif strcmp(pars.method,'PUM')
            l = legend('Patches','Node points','Eval points','Interpreter','latex','location','south','Orientation','Horizontal');
            l.NumColumns = 2;
        end
        l.FontSize = 20; 
    else
        if strcmp(pars.method,'FD')
            l = legend('Node points','Interpreter','latex','location','south','Orientation','Horizontal');
        elseif strcmp(pars.method,'PUM')
            l = legend('Patches','Node points','Interpreter','latex','location','south','Orientation','Horizontal');
        end
        l.FontSize = 20; 
    end

    % Plotting 
    plotSolution(ue,uExact,xe,dataY.bnd,pars.dim,pars.geom);
    %
    % Operator and solution l2 errors
    %
    lapNumeric = Lglobal(dataY.inner,:)*ucExact;
    evalNumeric = Eglobal(dataY.bnd,:)*ucExact;
    laplaceError = norm(lapAnalytic-lapNumeric,2)/norm(lapAnalytic,2);
    bndError = norm(evalNumeric-bndAnalytic,2)/(norm(bndAnalytic,2) + double(max(abs(bndAnalytic))==0));
    disp(['PDE error = ', num2str(l2Error)]);
    disp(['Boundary error = ', num2str(bndError)]);
    disp(['Laplace error = ', num2str(laplaceError)]);
end
end

