% -------------------------------------------------------------------------
% RBFsolver.m -- Numerically compute the solution to a Poisson problem using 
%                RBF based methods. Modules cpp, np, sp required (See
%                importModules) if pars.cppOn = 1.
% Inputs         -- pars                -> Parameter structure, see setPars.
%                   cpp                 -> Shared python library generated using
%                                          pybind11. See importModules.
%                   np                  -> Numpy library. See importModules.
%                   sp                  -> Scipy library. See importModules.
% Outputs        -- results             -> Output structure with fields:
%                   results.u           -> Double Nxd array. Numerical solution on centre (xc) points.
%                   results.ue          -> Double Mxd array. Numerical solution on eval (xe) points.
%                   results.ucExact     -> Double Nxd array. True solution on centre (xc) points.
%                   results.uExact      -> Double Mxd array. True solution on eval (xe) points.
%                   results.l2Error     -> Double. Discrete numerical error in the l2 norm.
%                   results.h           -> Double. Fill distance between node points.
%                   results.solveTime   -> Double. Solve time in seconds.
% Syntax         -- pars = setPars;
%                   [cpp, np, sp] = importModules;
%                   results = RBFsolver(pars, cpp, np, sp)   
%
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function results = RBFsolver(pars,varargin)
%
% Check input
%
if nargin > 4 || nargin < 1
    error("RBFsolver:IncorrectType","RBFsolver: Incorrect number of inputs.")
elseif ~isPars(pars)
    error("RBFsolver:IncorrectType","RBFsolver: pars should be a structure as given by setPars function.")
end
if pars.cppOn == 1 
    if nargin == 4 && (strcmp(underlyingType(varargin{1}),'py.module') && strcmp(underlyingType(varargin{2}),'py.module') && strcmp(underlyingType(varargin{3}),'py.module')) 
        [cpp, np, sp] = varargin{1:3};
    else
        error("RBFsolver:IncorrectType","RBFsolver: If parameter cppOn = 1, function requires python modules cpp, np, sp. Run [cpp, np, sp] = importModules and then RBFsolver(pars,cpp,np,sp).");
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
        pars.M = pars.N;                                      % Number of evaluation points (Y)
        dataY = dataX;
        xe = dataY.nodes;
    else
        pars.M = ceil(pars.q*pars.N);                                    % Number of evaluation points (Y)
        [dataY] = getPts(pars.geom,pars.M,n,C,R,"fitted",0);
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
        pars.N = size(xc,1);
        dataX.nodes = xc;
    end
    ptch = knnsearch(xc,xe,'K',1); % This is a list of closest centre points to eval points (centres)
    pars.n = n;
elseif strcmp(pars.method,'PUM')
    %
    % Get patches, evaluation (dataY) and centre (dataX) points
    %
    [ptch, dataX, dataY] = getPtch(C,R,n,pars);
    xc = dataX.nodes;
    xe = dataY.nodes;
    pars.N = size(xc,1);
    pars.M = size(xe,1);
    pars.n = n;
end
%
% Constructing global approximation to evaluation and Laplace operators
%
[L,B, Eglobal, Lglobal] = conGlobMat(pars,dataY,dataX,ptch);
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
    B = zeros(0,pars.N);
end
%
% Solution on centre points, evaluated on evaluation point set
%
A = [L; B];
if pars.cppOn   % Use c++ implementation
    checkMemory(A,F,pars.memTol);
    [BI, BJ, BV] = find(A);
    row  = np.array(BI - 1).reshape(int32(-1));
    col  = np.array(BJ - 1).reshape(int32(-1));
    data = np.array(BV,pyargs('dtype', 'float64')).reshape(int32(-1));
    pyB = sp.csc_matrix({data,{row, col}}, {int32(size(A,1)), int32(size(A,2))});
    u = double(cpp.sparseSolve(pyB,np.array(F,pyargs('dtype', 'float64')),pars.debug))';
    if pars.debug    % Save time to solve for comparisons
        sTime = u(end);
        u(end) = [];
        if pars.display
            disp(['RBFsolver: c++ sparse QR solver solve time was ', num2str(sTime), 's']);
        end
    end
else
    if pars.debug    % Save time to solve for comparisons
        tic
        u = A\F;
        sTime = toc;
        if pars.display
            disp(['RBFsolver: Matlab sparse QR solver solve time was ', num2str(sTime), 's']);
        end
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
% Plot points (evaluation and centre)
%
if pars.display
    if strcmp(pars.method, 'PUM')
        plotPtch(ptch,pars.geom,C,R)
        hold on
    end
    if pars.dim == 1
        evalPtPlot = plot(xe,zeros(size(xe,1),1),'b.'); hold on;
        centerPtPlot = plot(xc,zeros(size(xc,1),1),'rx');
        ax = gca;
        ax.LineWidth = 1.5;
        ax.FontSize = 18;
        xlabel("x",'Interpreter','latex','FontSize',20)
        xlim([-1.5 1.5]);
    elseif pars.dim == 2
        centerPtPlot = plot(xc(:,1),xc(:,2),'rx','MarkerSize',8,'LineWidth',3); hold on;
        if ~strcmp(pars.mode,"collocation")
            evalPtPlot = plot(xe(:,1),xe(:,2),'b.','MarkerSize',8);
        end
        ax = gca;
        ax.LineWidth = 1.5;
        ax.FontSize = 18;
        xlabel("x",'Interpreter','latex','FontSize',20)
        ylabel("y",'Interpreter','latex','FontSize',20,'Rotation',0)
        axis equal
        xlim([-1.5 1.5]);
        ylim([-1.5 1.5]);
    elseif pars.dim == 3
        evalPtPlot = plot3(xe(:,1),xe(:,2),xe(:,3),'b.'); hold on;
        centerPtPlot = plot3(xc(:,1),xc(:,2),xc(:,3),'rx');
        ax = gca;
        ax.LineWidth = 1.5;
        ax.FontSize = 18;
        xlabel("x",'Interpreter','latex','FontSize',20)
        ylabel("y",'Interpreter','latex','FontSize',20)
        ylabel("z",'Interpreter','latex','FontSize',20,'Rotation',0)
        axis equal
        xlim([-1.5 1.5]);
        ylim([-1.5 1.5]);
        zlim([-1.5 1.5]);
    end
    %
    % Plotting solution
    %
    plotSolution(ue,uExact,xe,dataY.bnd,pars.dim,pars.geom);
    %
    % Operator and solution l2 errors
    %
    lapNumeric = Lglobal(dataY.inner,:)*ucExact;
    evalNumeric = Eglobal(dataY.bnd,:)*ucExact;
    laplaceError = norm(lapAnalytic-lapNumeric,2)/norm(lapAnalytic,2);
    bndError = norm(evalNumeric-bndAnalytic,2)/(norm(bndAnalytic,2) + double(max(abs(bndAnalytic))==0));
    disp(['RBFsolver: PDE error = ', num2str(l2Error)]);
    disp(['RBFsolver: Boundary error = ', num2str(bndError)]);
    disp(['RBFsolver: Laplace error = ', num2str(laplaceError)]);
end
%
% Export results
%
results.u = u;
results.ue = ue;
results.ucExact = ucExact;
results.uExact = uExact;
results.l2Error = l2Error;
results.h = h;
if pars.debug
    results.solveTime = sTime;
end
end

