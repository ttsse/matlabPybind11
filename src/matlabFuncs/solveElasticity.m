% ----------------------------------------------------------------------------------
% solveElasticity.m -- form the derivative approximations based on the
%                      discretisation in Patches, formulate the sparse 
%                      least squares system and solve it using sparse QR.
% Copyright (c) 2023 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% ----------------------------------------------------------------------------------

function [pars, Result] = solveElasticity(geom,X_points,Y_points,Patches,pars)


% Extract Parameters

xb = Y_points.xb;
xi = Y_points.xi;
nb = Y_points.nb;
plist_i = Y_points.plist_i;
plist_b = Y_points.plist_b;
hy_b = Y_points.hy_b;
hy_i = Y_points.hy_i;
dim = size(xb,2);

xc = X_points.xc;
h = X_points.h;

Ni = size(xi,1);
Nb = size(xb,1);

C = Patches.C;
R = Patches.R;
Z = Patches.Z;
T = Patches.T;
psi = Patches.psi;
PDEapprox = Patches.PDEapprox;

clear X_points Y_points Patches 

% Lame parameters
lambda = (pars.E*pars.nu)/((1+pars.nu)*(1-2*pars.nu)); 
mu = pars.E/(2*(1+pars.nu));

mode = pars.mode;

% Select boundary conditions
switch mode
    case 'manufactured'
        cell_sz = length(pars.u_exact);
        bp_i = [min(xi(:,1)) pars.break_points max(xi(:,1))+1e-6];
        bp_b = [min(xb(:,1)) pars.break_points max(xb(:,1))+1e-6];
    case 'physical'
        pressure_b = pars.pressure_bottom; % bottom pressure
        pressure_t = pars.pressure_top; % top pressure
        body_force.x = pars.body_force.x;
        body_force.y = pars.body_force.y;
        body_force.z = pars.body_force.z;
end


%
% Construct discretized Elasticity matrices
% 

% Compute weight matrices
pu_i = cylPUWeights(psi,C,R,Z,T,xi,plist_i);
pu_b = cylPUWeights(psi,C,R,Z,T,xb,plist_b);
op = 'H'; % compute hessian
Hi=computeGlobalOp(PDEapprox,C,R,Z,T,op,xi,plist_i,pu_i);
op = 'J'; % compute jacobian
Jb=computeGlobalOp(PDEapprox,C,R,Z,T,op,xb,plist_b,pu_b);
op = '0'; % compute interpolation operators
Eb=computeGlobalOp(PDEapprox,C,R,Z,T,op,xb,plist_b,pu_b);
Ei=computeGlobalOp(PDEapprox,C,R,Z,T,op,xi,plist_i,pu_i);

%
% Compute differrors for each derivative
% 
differror = computeDifferror(xc,xi,xb,Ei,Eb,Jb,Hi);

% Elasticity matrix (inner points) 
ELi = [(lambda + 2*mu)*Hi{1,1} + mu*(Hi{2,2}+Hi{3,3}) (lambda + mu)*Hi{1,2} (lambda + mu)*Hi{1,3}; ...
       (lambda + mu)*Hi{2,1} (lambda + 2*mu)*Hi{2,2} + mu*(Hi{1,1} + Hi{3,3}) (lambda + mu)*Hi{2,3}; ...
       (lambda + mu)*Hi{3,1} (lambda + mu)*Hi{3,2} (lambda + 2*mu)*Hi{3,3} + mu*(Hi{1,1} + Hi{2,2})]; 

% Dirichlet BC
ELd = [Eb spalloc(size(Eb,1),size(Eb,2),0) spalloc(size(Eb,1),size(Eb,2),0); ...
       spalloc(size(Eb,1),size(Eb,2),0) Eb spalloc(size(Eb,1),size(Eb,2),0); ...
       spalloc(size(Eb,1),size(Eb,2),0) spalloc(size(Eb,1),size(Eb,2),0) Eb];    

%
% Note indices of boundary points to apply Boundary Conditions
%


    switch mode
        case 'manufactured'
            syms x1 x2 x3
            X = [x1 x2 x3];

            for i = 1:cell_sz
                u_exact{i} = pars.u_exact{i}';
                epsilon{i} = ([diff(u_exact{i}(1),x1), 0.5*(diff(u_exact{i}(2),x1) + diff(u_exact{i}(1),x2)), 0.5*(diff(u_exact{i}(3),x1) + diff(u_exact{i}(1),x3)); ...
                        0.5*(diff(u_exact{i}(1),x2) + diff(u_exact{i}(2),x1)), diff(u_exact{i}(2),x2), 0.5*(diff(u_exact{i}(3),x2) + diff(u_exact{i}(2),x3)); ...
                        0.5*(diff(u_exact{i}(1),x3) + diff(u_exact{i}(3),x1)), 0.5*(diff(u_exact{i}(3),x2) + diff(u_exact{i}(2),x3)), diff(u_exact{i}(3),x3);]);
                sigma{i} = lambda.*trace(epsilon{i}).*eye(3,3) + 2.*mu.*epsilon{i};
                f_exact{i} = -[divergence(sigma{i}(:,1),X); divergence(sigma{i}(:,2),X); divergence(sigma{i}(:,3),X)];
            end
            
            % Dirichlet
            d1_xyz = find(xb(:,1) == -(geom.l/2));
            
            d2_xyz = find(xb(:,1) == (geom.l/2));
            
            % Get Neumann points
            N1_xyz = find(xb(:,3) == -(geom.h/2) & xb(:,1) ~= -(geom.l/2) & xb(:,1) ~= (geom.l/2));
            N_n1 = repmat([0, 0, -1],length(N1_xyz),1);

            N2_xyz = find(xb(:,3) == (geom.h/2) & xb(:,1) ~= -(geom.l/2) & xb(:,1) ~= (geom.l/2));
            N_n2 = repmat([0, 0, 1],length(N2_xyz),1);

            N3_xyz = find(xb(:,2) == -(geom.w/2) & xb(:,1) ~= -(geom.l/2) & xb(:,1) ~= (geom.l/2) & xb(:,3) ~= (geom.h/2) & xb(:,3) ~= -(geom.h/2));
            N_n3 = repmat([0, -1, 0],length(N3_xyz),1);

            N4_xyz = find(xb(:,2) == (geom.w/2) & xb(:,1) ~= -(geom.l/2) & xb(:,1) ~= (geom.l/2) & xb(:,3) ~= (geom.h/2) & xb(:,3) ~= -(geom.h/2));
            N_n4 = repmat([0, 1, 0],length(N4_xyz),1);
                            
        case 'physical'
            d1_xyz = find(xb(:,1) == -(geom.l/2));

            d2_xyz = find(xb(:,1) == (geom.l/2));

            N1_xyz = find(xb(:,3) == -(geom.h/2) & xb(:,1) ~= -(geom.l/2) & xb(:,1) ~= (geom.l/2));
            N_n1 = repmat([0, 0, -1],size(N1_xyz,1),1);

            N2_xyz = find(xb(:,3) == (geom.h/2) & xb(:,1) ~= -(geom.l/2) & xb(:,1) ~= (geom.l/2));
            N_n2 = repmat([0, 0, 1],size(N2_xyz,1),1);

            N3_xyz = find(xb(:,2) == -(geom.w/2) & xb(:,1) ~= -(geom.l/2) & xb(:,1) ~= (geom.l/2) & xb(:,3) ~= (geom.h/2) & xb(:,3) ~= -(geom.h/2));
            N_n3 = repmat([0, -1, 0],size(N3_xyz,1),1);

            N4_xyz = find(xb(:,2) == (geom.w/2) & xb(:,1) ~= -(geom.l/2) & xb(:,1) ~= (geom.l/2) & xb(:,3) ~= (geom.h/2) & xb(:,3) ~= -(geom.h/2));
            N_n4 = repmat([0, 1, 0],size(N4_xyz,1),1);

            p_top = pressure_t.*N_n2;
    end

%-----------------------------------------------------------
% Compute the right hand side including scaling and scale the 
% different blocks of the PDE operator matrix in the same way
%-----------------------------------------------------------

switch mode
    case 'manufactured'
        ELn_func = @(n1,n2,n3,i) [(lambda + 2*mu)*n1*Jb{1}(i,:) + mu*(n2*Jb{2}(i,:)+n3*Jb{3}(i,:)) lambda*n1*Jb{2}(i,:) + mu*n2*Jb{1}(i,:) lambda*n1*Jb{3}(i,:) + mu*n3*Jb{1}(i,:); ...
                                lambda*n2*Jb{1}(i,:) + mu*n1*Jb{2}(i,:) (lambda + 2*mu)*n2*Jb{2}(i,:) + mu*(n1*Jb{1}(i,:)+n3*Jb{3}(i,:)) lambda*n2*Jb{3}(i,:) + mu*n3*Jb{2}(i,:); ...
                                lambda*n3*Jb{1}(i,:) + mu*n1*Jb{3}(i,:) lambda*n3*Jb{2}(i,:) + mu*n2*Jb{3}(i,:) (lambda + 2*mu)*n3*Jb{3}(i,:) + mu*(n1*Jb{1}(i,:)+n2*Jb{2}(i,:))]; 
        ELb = [];
        ELb = [ELb; ELn_func(spdiags(N_n1(:,1),0,length(N_n1),length(N_n1)),...
            spdiags(N_n1(:,2),0,length(N_n1),length(N_n1)),...
            spdiags(N_n1(:,3),0,length(N_n1),length(N_n1)),N1_xyz)...
            .*pars.scaling.neumann(h,hy_b,dim)];
        ELb = [ELb; ELn_func(spdiags(N_n2(:,1),0,length(N_n2),length(N_n2)),...
            spdiags(N_n2(:,2),0,length(N_n2),length(N_n2)),...
            spdiags(N_n2(:,3),0,length(N_n2),length(N_n2)),N2_xyz)...
            .*pars.scaling.neumann(h,hy_b,dim)];    
        ELb = [ELb; ELn_func(spdiags(N_n3(:,1),0,length(N_n3),length(N_n3)),...
            spdiags(N_n3(:,2),0,length(N_n3),length(N_n3)),...
            spdiags(N_n3(:,3),0,length(N_n3),length(N_n3)),N3_xyz)...
            .*pars.scaling.neumann(h,hy_b,dim)];  
        ELb = [ELb; ELn_func(spdiags(N_n4(:,1),0,length(N_n4),length(N_n4)),...
            spdiags(N_n4(:,2),0,length(N_n4),length(N_n4)),...
            spdiags(N_n4(:,3),0,length(N_n4),length(N_n4)),N4_xyz)...
            .*pars.scaling.neumann(h,hy_b,dim)];  
        ELb = [ELb; ELd([d1_xyz; d1_xyz + Nb; d1_xyz + 2*Nb],:).*pars.scaling.dirichlet(h,hy_b,dim)];
        ELb = [ELb; ELd([d2_xyz; d2_xyz + Nb; d2_xyz + 2*Nb],:).*pars.scaling.dirichlet(h,hy_b,dim)];
    case 'physical'
                ELn_func = @(n1,n2,n3,i) [(lambda + 2*mu)*n1*Jb{1}(i,:) + mu*(n2*Jb{2}(i,:)+n3*Jb{3}(i,:)) lambda*n1*Jb{2}(i,:) + mu*n2*Jb{1}(i,:) lambda*n1*Jb{3}(i,:) + mu*n3*Jb{1}(i,:); ...
                                        lambda*n2*Jb{1}(i,:) + mu*n1*Jb{2}(i,:) (lambda + 2*mu)*n2*Jb{2}(i,:) + mu*(n1*Jb{1}(i,:)+n3*Jb{3}(i,:)) lambda*n2*Jb{3}(i,:) + mu*n3*Jb{2}(i,:); ...
                                        lambda*n3*Jb{1}(i,:) + mu*n1*Jb{3}(i,:) lambda*n3*Jb{2}(i,:) + mu*n2*Jb{3}(i,:) (lambda + 2*mu)*n3*Jb{3}(i,:) + mu*(n1*Jb{1}(i,:)+n2*Jb{2}(i,:))]; 
                ELb = [];
                ELb = [ELb; ELn_func(spdiags(N_n1(:,1),0,size(N_n1,1),size(N_n1,1)),...
                    spdiags(N_n1(:,2),0,size(N_n1,1),size(N_n1,1)),...
                    spdiags(N_n1(:,3),0,size(N_n1,1),size(N_n1,1)),N1_xyz)...
                    .*pars.scaling.neumann(h,hy_b,dim)];
                ELb = [ELb; ELn_func(spdiags(N_n2(:,1),0,size(N_n2,1),size(N_n2,1)),...
                    spdiags(N_n2(:,2),0,size(N_n2,1),size(N_n2,1)),...
                    spdiags(N_n2(:,3),0,size(N_n2,1),size(N_n2,1)),N2_xyz)...
                    .*pars.scaling.neumann(h,hy_b,dim)];
                ELb = [ELb; ELn_func(spdiags(N_n3(:,1),0,size(N_n3,1),size(N_n3,1)),...
                    spdiags(N_n3(:,2),0,size(N_n3,1),size(N_n3,1)),...
                    spdiags(N_n3(:,3),0,size(N_n3,1),size(N_n3,1)),N3_xyz)...
                    .*pars.scaling.neumann(h,hy_b,dim)];  
                ELb = [ELb; ELn_func(spdiags(N_n4(:,1),0,size(N_n4,1),size(N_n4,1)),...
                    spdiags(N_n4(:,2),0,size(N_n4,1),size(N_n4,1)),...
                    spdiags(N_n4(:,3),0,size(N_n4,1),size(N_n4,1)),N4_xyz)...
                    .*pars.scaling.neumann(h,hy_b,dim)]; 
                ELb = [ELb; ELd([d1_xyz; d1_xyz + Nb; d1_xyz + 2*Nb],:).*pars.scaling.dirichlet(h,hy_b,dim)];
                ELb = [ELb; ELd([d2_xyz; d2_xyz + Nb; d2_xyz + 2*Nb],:).*pars.scaling.dirichlet(h,hy_b,dim)];
end

ELb = sparse(ELb);
ELi = ELi.*pars.scaling.interior(h,hy_i,dim);  
K = [-ELi; ELb]; % LHS

f_i = zeros(Ni*3,1);
f_b = zeros(Nb*3,1);


switch mode
    case 'manufactured'
            clear f_b
            f_i = p_func(f_exact,xi,bp_i,cell_sz).*pars.scaling.interior(h,hy_i,dim);
            traction = pagemtimes(p_func(sigma,xb(N1_xyz,:),bp_b,cell_sz),reshape(N_n1',[3 1 length(N_n1)]));
            f_b_N1 = reshape(permute(traction,[3 1 2]),[],1).*pars.scaling.neumann(h,hy_b,dim);
            traction = pagemtimes(p_func(sigma,xb(N2_xyz,:),bp_b,cell_sz),reshape(N_n2',[3 1 length(N_n2)]));
            f_b_N2 = reshape(permute(traction,[3 1 2]),[],1).*pars.scaling.neumann(h,hy_b,dim);
            traction = pagemtimes(p_func(sigma,xb(N3_xyz,:),bp_b,cell_sz),reshape(N_n3',[3 1 length(N_n3)]));
            f_b_N3 = reshape(permute(traction,[3 1 2]),[],1).*pars.scaling.neumann(h,hy_b,dim);
            traction = pagemtimes(p_func(sigma,xb(N4_xyz,:),bp_b,cell_sz),reshape(N_n4',[3 1 length(N_n4)]));
            f_b_N4 = reshape(permute(traction,[3 1 2]),[],1).*pars.scaling.neumann(h,hy_b,dim);
            
            f_b_d1 = p_func(u_exact,xb(d1_xyz,:),bp_b,cell_sz).*pars.scaling.dirichlet(h,hy_b,dim);
            f_b_d2 = p_func(u_exact,xb(d2_xyz,:),bp_b,cell_sz).*pars.scaling.dirichlet(h,hy_b,dim);
            f_b = [f_b_N1; f_b_N2; f_b_N3; f_b_N4; f_b_d1; f_b_d2];
    case 'physical'
        clear f_b
            f_i(1:Ni,1) = body_force.x.*pars.scaling.interior(h,hy_i,dim);
            f_i(Ni+1:2*Ni,1) = body_force.y.*pars.scaling.interior(h,hy_i,dim);
            f_i(2*Ni+1:3*Ni,1) = body_force.z.*pars.scaling.interior(h,hy_i,dim);                                
            f_b_N1 = zeros(3*size(N1_xyz,1),1).*pars.scaling.neumann(h,hy_b,dim);
            f_b_N2 = [p_top(:,1); p_top(:,2); p_top(:,3)].*pars.scaling.neumann(h,hy_b,dim);
            f_b_N3 = zeros(3*size(N3_xyz,1),1).*pars.scaling.neumann(h,hy_b,dim);
            f_b_N4 = zeros(3*size(N4_xyz,1),1).*pars.scaling.neumann(h,hy_b,dim);
            f_b_d1 = zeros(3*size(d1_xyz,1),1).*pars.scaling.dirichlet(h,hy_b,dim);
            f_b_d2 = zeros(3*size(d2_xyz,1),1).*pars.scaling.dirichlet(h,hy_b,dim);
            f_b = [f_b_N1; f_b_N2; f_b_N3; f_b_N4; f_b_d1; f_b_d2];
end


f = [f_i; f_b]; % RHS

%-----------------------------------------------------------
% Solve the problem 
%-----------------------------------------------------------

% Clear extra matrices to free space
clear f_i f_b ELb ELd ELi ELn Hi n1 n2 n3 nb ELb_d ELb_n f_b_N1 f_b_N2 f_b_N3 f_b_N4


[KI, KJ, KV] = find(K);
pyK = py.scipy.sparse.csc_matrix({KV, {uint64(KI-1) uint64(KJ-1)}}, {uint64(size(K,1)), uint64(size(K,2))});

rbfsol = double(py.cppFuncs.sparseSolve(pyK,py.numpy.array(f,copy=false,order='F')))';

residual = norm((K*rbfsol) - f,2);

tic
rbfsol2 = K\f;
toc

if ~setdiff(rbfsol,rbfsol2)
    error("Incorrect solution")
end

clear K f

%---------------------------------------------
% Evaluate displacement and stress on Y points
%---------------------------------------------

xall = [xi; xb];
switch mode
    case 'manufactured'
        u_actual = p_func(u_exact,xall,bp_b,cell_sz);
        evalu = [[Ei; Eb]*rbfsol(1:size(Ei,2)); ...
         [Ei; Eb]*rbfsol(size(Ei,2)+1:2*size(Ei,2)); ...
         [Ei; Eb]*rbfsol(2*size(Ei,2)+1:3*size(Ei,2))];
        L2_error = norm(u_actual - evalu)/norm(u_actual);
        Linf_error = norm(u_actual - evalu,Inf)/norm(u_actual,Inf);
        u_actual = [u_actual(1:Nb+Ni) u_actual(Nb+Ni+1:2*(Nb+Ni)) u_actual(2*(Nb+Ni)+1:3*(Nb+Ni))];
        Result.u_exact = u_actual;
    case 'physical'
        evalu = [[Ei; Eb]*rbfsol(1:size(Ei,2)); ...
         [Ei; Eb]*rbfsol(size(Ei,2)+1:2*size(Ei,2)); ...
         [Ei; Eb]*rbfsol(2*size(Ei,2)+1:3*size(Ei,2))];
end

evalu = [evalu(1:Nb+Ni) evalu(Nb+Ni+1:2*(Nb+Ni)) evalu(2*(Nb+Ni)+1:3*(Nb+Ni))];

x_new = xall + evalu;
xi_n = x_new(1:length(xi),:);
xb_n = x_new(length(xi)+1:length(xi)+length(xb),:);

% Plot solution (movement)
figure()
plot3(xi_n(:,1),xi_n(:,2),xi_n(:,3),'r.'), hold on
plot3(xb_n(:,1),xb_n(:,2),xb_n(:,3),'r.')
hold on
plot3(xi(:,1),xi(:,2),xi(:,3),'g.'), hold on
plot3(xb(:,1),xb(:,2),xb(:,3),'g.')
axis equal

%----------------------------------------------------------
% Export Results
%----------------------------------------------------------
Result.xb = xb; % boundary y points
Result.xi = xi; % inner y points
Result.xall = xall; % y points
Result.xb_n = xb_n; % boundary y points solution
Result.xi_n = xi_n; % inner y points solution
Result.sol_y = evalu; % solution on y points
Result.sol_x = rbfsol; % solution on x points
Result.residual = residual;
Result.differror = differror;
Result.h = h;

switch mode 
    case 'manufactured'
        Result.L2 = L2_error;
        Result.Linf = Linf_error; 
end
end