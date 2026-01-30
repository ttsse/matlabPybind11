%
% Plotting routines
%
function [] = plotSolution(uNumeric,uAnalytic,x,idB,dim,geom)
    err = uNumeric-uAnalytic;
    if dim == 1
        [x,id] = sort(x);
        figure()
        plot(x,uAnalytic(id),'r-'); 
        hold on;
        plot(x,uNumeric(id),'bo');
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("$$u$$","Interpreter","latex","FontSize",24,'Rotation',0)
        legend('$u_E$','$u_N$','FontSize',18,'Interpreter','latex');
        grid on
    
        figure()
        plot(x,err(id),'b-.');
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("$e_{rel}$","Interpreter","latex","FontSize",24,'Rotation',0)
        grid on
    elseif dim == 2
        T = delaunay(x(:,1),x(:,2));

        figure()
        G=trisurf(T,x(:,1),x(:,2),uAnalytic);
        hold on
        if strcmp(geom,"ball")
            plot(x([idB; idB(1)],1),x([idB; idB(1)],2),'k-',"LineWidth",1.5)
        elseif strcmp(geom,"cube")
            plot([max(x(idB,1)) max(x(idB,1)) min(x(idB,1)) min(x(idB,1)) max(x(idB,1))],...
                 [max(x(idB,2)) min(x(idB,2)) min(x(idB,2)) max(x(idB,2)) max(x(idB,2))],'k-',"LineWidth",1.5)
        end
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("$$u_E$$","Interpreter","latex","FontSize",24,'Rotation',0)
        set(G,'EdgeColor','none')
        shading interp
    
        figure()
        G=trisurf(T,x(:,1),x(:,2),uNumeric);
        hold on
        if strcmp(geom,"ball")
            plot(x([idB; idB(1)],1),x([idB; idB(1)],2),'k-',"LineWidth",1.5)
        elseif strcmp(geom,"cube")
            plot([max(x(idB,1)) max(x(idB,1)) min(x(idB,1)) min(x(idB,1)) max(x(idB,1))],...
                 [max(x(idB,2)) min(x(idB,2)) min(x(idB,2)) max(x(idB,2)) max(x(idB,2))],'k-',"LineWidth",1.5)
        end
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("$$u$$","Interpreter","latex","FontSize",24,'Rotation',0)
        set(G,'EdgeColor','none')
        shading interp
    
        figure()
        G=trisurf(T,x(:,1),x(:,2),err);
        hold on
        if strcmp(geom,"ball")
            plot(x([idB; idB(1)],1),x([idB; idB(1)],2),'k-',"LineWidth",1.5)
        elseif strcmp(geom,"cube")
            plot([max(x(idB,1)) max(x(idB,1)) min(x(idB,1)) min(x(idB,1)) max(x(idB,1))],...
                 [max(x(idB,2)) min(x(idB,2)) min(x(idB,2)) max(x(idB,2)) max(x(idB,2))],'k-',"LineWidth",1.5)
        end
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("$$u-u_E$$","Interpreter","latex","FontSize",24)
        set(G,'EdgeColor','none')
        shading interp
    elseif dim == 3 && ~strcmp(geom,"cube")
        
        bndPtCloud = pointCloud(x(idB,:));
        surfMesh = pc2surfacemesh(bndPtCloud,"ball-pivot");
        x = surfMesh.Vertices;
        T = surfMesh.Faces;
        
        figure()
        G=trisurf(T,x(:,1),x(:,2),x(:,3),uAnalytic(idB));
        axis equal
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("z","Interpreter","latex","FontSize",24,"Rotation",0)
        set(G,'EdgeColor','none')
        cb = colorbar;
        ylabel(cb,"$$u_E$$","Interpreter","latex","FontSize",24,'Rotation',0)
        shading interp
        
        figure()
        G=trisurf(T,x(:,1),x(:,2),x(:,3),uNumeric(idB));
        axis equal
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("z","Interpreter","latex","FontSize",24,"Rotation",0)
        set(G,'EdgeColor','none')
        cb = colorbar;
        ylabel(cb,"$$u$$","Interpreter","latex","FontSize",24,'Rotation',0)
        shading interp
    
        figure()
        G=trisurf(T,x(:,1),x(:,2),x(:,3),err(idB));
        axis equal
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("z","Interpreter","latex","FontSize",24,"Rotation",0)
        set(G,'EdgeColor','none')
        cb = colorbar;
        ylabel(cb,"$$u-u_E$$","Interpreter","latex","FontSize",24,'Rotation',90)
        shading interp
    else
        disp("No plotting availablein 3D for this geometry")
    end
end