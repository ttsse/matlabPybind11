% Function that puts Halton points in a cuboid with size lxwxh with mean
% distance between points, h_p

function [bdata geom] = makeFEMRectangle(l,w,h,display)
        
    filename = 'cuboidGeom';
               
    A = readmatrix(filename);

    for i = 1:size(A,1)/3
        x(i,:) = A(3*i,1:3);
    end
    
    H_top = x(find(x(:,3) == h/2),:);  n_top = [zeros(size(H_top,1),1), zeros(size(H_top,1),1), ones(size(H_top,1),1)];
    H_bot = x(find(x(:,3) == -h/2),:); n_bot = [zeros(size(H_top,1),1), zeros(size(H_top,1),1), -ones(size(H_top,1),1)];
    H_s1 = x(find(x(:,1) == l/2 & x(:,2) ~= w/2 & x(:,3) ~= h/2 & x(:,2) ~= -w/2 & x(:,3) ~= -h/2),:); n_s1 = [ones(size(H_s1,1),1), zeros(size(H_s1,1),1), zeros(size(H_s1,1),1)];
    H_s2 = x(find(x(:,1) == -l/2 & x(:,2) ~= w/2 & x(:,3) ~= h/2 & x(:,2) ~= -w/2 & x(:,3) ~= -h/2),:); n_s2 = [-ones(size(H_s2,1),1), zeros(size(H_s2,1),1), zeros(size(H_s2,1),1)];
    H_s3 = x(find(x(:,2) == w/2 & x(:,1) ~= l/2 & x(:,3) ~= h/2 & x(:,1) ~= -l/2 & x(:,3) ~= -h/2),:); n_s3 = [zeros(size(H_s3,1),1), ones(size(H_s3,1),1), zeros(size(H_s3,1),1)];
    H_s4 = x(find(x(:,2) == -w/2 & x(:,1) ~= l/2 & x(:,3) ~= h/2 & x(:,1) ~= -l/2 & x(:,3) ~= -h/2),:); n_s4 = [zeros(size(H_s4,1),1), -ones(size(H_s4,1),1), zeros(size(H_s4,1),1)];
    H_b = [H_top; H_bot; H_s1; H_s2; H_s3; H_s4];
    n = [n_top; n_bot; n_s1; n_s2; n_s3; n_s4];
    
    H_inside = setdiff(x,H_b,'rows');
    
    H = [H_inside; H_b];
    
    bdata.nodes = H_b;
    bdata.normal = n;
    bdata.outer = find(bdata.nodes(:,3) == -h/2);
    bdata.inner = find(bdata.nodes(:,3) ~= -h/2);
    bdata.xi = H_inside;
    bdata.xb_i = bdata.nodes(bdata.inner,:);
    bdata.xb_o = bdata.nodes(bdata.outer,:);
    geom.l = l;
    geom.w = w;
    geom.h = h;
    geom.psi = 'bump';
    
    % Plotting
    if (display)
        figure();
        plot3(H(:,1),H(:,2),H(:,3),'.');
        axis equal

        figure();
        quiver3(bdata.nodes(:,1),bdata.nodes(:,2),bdata.nodes(:,3),n(:,1),n(:,2),n(:,3));
        axis equal
    end

end