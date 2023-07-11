% -------------------------------------------------------------------------
% RBFsolver.m -- Numerically compute the solution to a Poisson problem using 
%                RBF based methods.
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function [l2Error, h] = RBFsolver(pars)

method = pars.method;   % FD or PUM
dim = pars.dim;                            % dim = 1,2 or 3
display = pars.display;                        % Plot solution
geom = pars.geom;                      % ball or cube
mode = pars.mode;                    % fitted, unfitted or collocation
bcMode = pars.bcMode;                    % strong or weak imposition of boundary conditions (only relevant for fitted)
scaling = pars.scaling;                        % Include scaling of the unfitted LS problem
mvCentres = pars.mvCentres;                      % Option to have a Y point on top of all X points inside the domain
q = pars.q;                              % Oversampling for LS methods
N = pars.N;                             % Number of centre points (X) in each patch - RBFPUM / Number of all centre points (X) - RBFFD 
P = pars.P;
ep = pars.ep;                           % Not relevant for 'r3' basis
phi = pars.phi;                      % Choice of basis 'r3', 'mq', 'gs', 'iq', 'rbfqr'
psi = pars.psi;
pdeg = pars.pdeg;                          % Polynomial extension, not relevant for 'rbfqr'
rbfdeg = pars.rbfdeg;
extCoeff = pars.extCoeff;
del = pars.del;
prob = pars.prob;

% 
% Point number for RBF-FD method choice
%
if strcmp(method,'FD')
    if pdeg == -1
        n = nchoosek(rbfdeg+dim-1,dim); % Stencil size in smooth RBF case
    else
        n = 2*nchoosek(pdeg+dim,dim);   % Stencil size in phs + poly case
    end
    
    if N < n && strcmp(method,'FD')
        warning("RBFFD: The number of centre points are less than the stencil points required for the given polynomial order. N is increased.")
        N = n;
    end
    extCoeff = extCoeff*(strcmp(mode,"unfitted"));
elseif strcmp(method,'PUM')
    n = nchoosek(pars.rbfdeg+pars.dim,pars.dim);
end

%
% Place N centre points and M evaluation points in geom with centre C and radius R
%
C = zeros(1,dim);
if strcmp(geom,'cube')
    R = 1*(dim)^(1/2);
elseif strcmp(geom,'ball')
    R = 1;
else
    error("RBFsolver: Requested geometry not implemented yet");
end

if strcmp(method,'FD')
    %
    % Get centre points (xc)
    %
    [dataX] = getPts(geom,N,n,C,R,mode,extCoeff);
    xc = dataX.nodes;                                     
    % 
    % Get evaluation points (xe)
    %
    if strcmp(mode,"collocation")       
        M = N;                                      % Number of evaluation points (Y)
        dataY = dataX;
        xe = dataY.nodes;
    else
        M = ceil(q*N);                                    % Number of evaluation points (Y)
        [dataY] = getPts(geom,M,n,C,R,"fitted",0);
        %
        % Move evaluation points inside (and on the boundary) to the closest center point.
        % Make sure not to move two evaluation points to the same center point. All
        % center points should be moved to.
        %
        if mvCentres && ~strcmp(mode,"collocation")
            dataY = movePts(dataX,dataY);
        end
        xe = dataY.nodes;
    end
    %
    % Ensure center points are not too far from boundary and make evaluation point - stencil list 
    %
    if strcmp(mode,"unfitted")
        xc = xc(unique(knnsearch(xc,xe,'K',ceil(extCoeff*n+eps(1)))),:);
        N = size(xc,1);
    end
    ptStencilList = knnsearch(xc,xe,'K',1);
elseif strcmp(method,'PUM')
    %
    % Get patches  
    %
    ptch = getPtch(geom,P,C,R,del);
    P = length(ptch.R);
    %
    % Get centre (xc)
    %
    if strcmp(mode,"unfitted")
        dataX.nodes = [];            % Centers
        for i = 1:P
            [ptch.xc(i)] = getPts("ball",n,0,ptch.C(i,:),ptch.R(i),"unfitted",0);
            xCglobalId{i} = [(i-1)*n+1:i*n]';
            dataX.nodes = [dataX.nodes; ptch.xc(i).nodes];
        end
        xc = dataX.nodes;
        [ptch.xc.globalId] = xCglobalId{:};
        if strcmp(geom,"ball")
            dataX.inner = find(sqrt(sum((xc-C).^2,2))<=R);
        elseif strcmp(geom,"cube")
            dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];  
            dataX.inner = find(all(abs(xc-C)<=dimLCoeff(dim)*R,2));
        end
        dataX.outer = setdiff(1:size(xc,1),dataX.inner);
        dataX.bnd = [];
    elseif strcmp(mode,"collocation")
        dataX = getPts(geom,n*P,0,C,R,"fitted",0);
        xc = dataX.nodes;
        for i = 1:P
            ptch.xc(i).globalId = find(sqrt(sum((xc - ptch.C(i,:)).^2,2)) <= ptch.R(i));
            ptch.xc(i).nodes = xc(ptch.xc(i).globalId,:);
        end
    elseif strcmp(mode,"fitted") 
        dataX = getPts(geom,n*P,0,C,R,"fitted",0);
        xc = dataX.nodes;
        for i = 1:P
            ptch.xc(i).globalId = find(sqrt(sum((xc - ptch.C(i,:)).^2,2)) <= ptch.R(i));
            ptch.xc(i).nodes = xc(ptch.xc(i).globalId,:);
        end
    end
    %
    % Get evaluation points (xe)
    %
    q = max(q*double(~strcmp(mode,"collocation")),1);
    M = ceil(n*P*q);  
    if ~strcmp(mode,"collocation")
        dataY = getPts(geom,M,0,C,R,"fitted",0);
        %
        % Move evaluation points inside (and on the boundary) to the closest center point.
        %
        if mvCentres && ~strcmp(mode,"collocation")
            dataY = movePts(dataX,dataY);
        end
        xe = dataY.nodes;
        for i = 1:P
            ptch.xe(i).globalId = find(sqrt(sum((xe - ptch.C(i,:)).^2,2)) <= ptch.R(i));
            ptch.xe(i).nodes = xe(ptch.xe(i).globalId,:);
        end
    else
        dataY = dataX;
        xe = xc;        % Centers
        ptch.xe = ptch.xc;
    end
    N = size(xc,1);
end
%
% Constructing global LS-RBF-FD approximation to evaluation and Laplace operators M x N
%
if strcmp(method,'FD')
    [L,B, Eglobal, Lglobal] = conGlobMat(dataY,M,N,n,prob,phi,psi,ep,pdeg,ptStencilList,xc,xe);
elseif strcmp(method,'PUM')
    [L,B, Eglobal, Lglobal] = conGlobMat(dataY,M,N,n,prob,phi,psi,ep,pdeg,ptch);
end
%
% Manufactured solution to construct forcing and BC
%
if dim == 1
    fun = @(x) sin(2.*pi.*x);
    lapFun = @(x) -4.*(pi.^2).*sin(2.*pi.*x);
    F = [lapFun(xe(dataY.inner,1)); fun(xe(dataY.bnd,1))];
    uExact = fun(xe(:,1));
    ucExact = fun(xc(:,1));
    lapAnalytic = lapFun(xe(dataY.inner,1));
    bndAnalytic = [0; 0];
elseif dim == 2
    fun = @(x,y) sin(2.*pi.*x.*y);
    lapFun = @(x,y) - 4.*(x.^2).*(pi.^2).*sin(2.*pi.*x.*y) - 4.*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y);
    F = [lapFun(xe(dataY.inner,1),xe(dataY.inner,2)); fun(xe(dataY.bnd,1),xe(dataY.bnd,2))];
    uExact = fun(xe(:,1),xe(:,2));
    ucExact = fun(xc(:,1),xc(:,2));
    lapAnalytic = lapFun(xe(dataY.inner,1),xe(dataY.inner,2));
    bndAnalytic = fun(xe(dataY.bnd,1),xe(dataY.bnd,2));
elseif dim == 3
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
if scaling && (strcmp(mode,"unfitted") || (strcmp(mode,"fitted") && strcmp(bcMode,"weak")))
    Lscale = sqrt(dataY.Vol/length(dataY.inner));
    Bscale = sqrt(dataY.Area/length(dataY.bnd))*(h.^-1.5);
    L = Lscale.*L; B = Bscale.*B;
    F(dataY.inner) = Lscale.*F(dataY.inner);
    F(dataY.bnd) = Bscale.*F(dataY.bnd);
end

if strcmp(mode,"fitted") && strcmp(bcMode,"strong")
    Fmod = Lglobal(:,dataX.bnd)*ucExact(dataX.bnd);
    F = F-Fmod;
    L = Lglobal(:,dataX.inner);
    B = zeros(0,N);
end
%
% Solution on centre points, evaluated on Y set
%
A = [L; B];
testMemory(nnz(A),M,pars.memTol);
tic
u1 = A\F;
toc
[BI, BJ, BV] = find(A);
pyB = py.scipy.sparse.csc_matrix({BV, {uint64(BI-1) uint64(BJ-1)}}, {uint64(size(A,1)), uint64(size(A,2))});
tic
u = double(py.cppFuncs.sparseSolve(pyB,py.numpy.array(F,copy=false,order='F')))';
toc
%
% Fix operators to compute error measures
%
if strcmp(mode,"fitted") && strcmp(bcMode,"strong")
    u = [u; ucExact(dataX.bnd)];
end
ue = Eglobal*u;
l2Error = norm(abs(ue-uExact),2)/norm(uExact,2);
%
% Displaying output 
%
if display
    if strcmp(method,'FD')
        figure()
        hold on;
    elseif strcmp(method, 'PUM')
        plotPtch(ptch,geom,C,R)
        hold on
    end
    if dim == 1
        evalPtPlot = plot(xe,zeros(size(xe,1),1),'b.');
        centerPtPlot = plot(xc,zeros(size(xc,1),1),'rx');
    elseif dim == 2
        centerPtPlot = plot(xc(:,1),xc(:,2),'rx','MarkerSize',8,'LineWidth',3);
        if ~strcmp(mode,"collocation")
            evalPtPlot = plot(xe(:,1),xe(:,2),'b.','MarkerSize',8);
        end
    elseif dim == 3
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
    if strcmp(mode,"unfitted")
        if strcmp(method,'FD')
            l = legend('Node points','Eval points','Interpreter','latex','location','south','Orientation','Horizontal');
        elseif strcmp(method,'PUM')
            l = legend('Patches','Node points','Eval points','Interpreter','latex','location','south','Orientation','Horizontal');
            l.NumColumns = 2;
        end
        l.FontSize = 20; 
    else
        if strcmp(method,'FD')
            l = legend('Node points','Interpreter','latex','location','south','Orientation','Horizontal');
        elseif strcmp(method,'PUM')
            l = legend('Patches','Node points','Interpreter','latex','location','south','Orientation','Horizontal');
        end
        l.FontSize = 20; 
    end

    % Plotting 
    plotSolution(ue,uExact,xe,dataY.bnd,dim,geom);
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

