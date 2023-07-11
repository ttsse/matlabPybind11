function [L,B, Eglobal, Lglobal] = conGlobMat(dataY,M,N,n,prob,phi,psi,ep,pdeg,ptStencilList,xc,xe)
    
    if ~strcmp(prob,'Poisson')
        error("conGlobMat: Requested problem is not currently implemented")
    end
    
    if ~isstruct(ptStencilList)
        relXc = unique(ptStencilList); % stencil centres that are used directly
        Eglobal = spalloc(M,N,M*n);
        Lglobal = spalloc(M,N,M*n);
        for i = 1:length(relXc)
            [idX,~] = knnsearch(xc,xc(relXc(i),:),'K',n);
            xcLoc = xc(idX,:); % stencil
            Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));
            idY = find(ptStencilList==relXc(i));
            E = RBFDiffMat(0,Psi,xe(idY,:));
            L = RBFDiffMat(1.5,Psi,xe(idY,:));
            
            Eglobal(idY,idX) = E + Eglobal(idY,idX);
            Lglobal(idY,idX) = L + Lglobal(idY,idX);
        end
        L = Lglobal(dataY.inner,:);
        B = Eglobal(dataY.bnd,:);
    else
        ptch = ptStencilList;
        dim = size(ptch.C,2);
        P = length(ptch.R);
        Eglobal = spalloc(M,N,M*n);
        Lglobal = spalloc(M,N,M*n);
        [w] = weights(psi,1.5,ptch);
        for i = 1:P
            Psi = RBFInterpMat(phi,pdeg,ep,ptch.xc(i).nodes,ptch.C(i,:),ptch.R(i));
            E = RBFDiffMat(0,Psi,ptch.xe(i).nodes);
            B = RBFDiffMat(1,Psi,ptch.xe(i).nodes);
            L = RBFDiffMat(1.5,Psi,ptch.xe(i).nodes);
            Eglobal(ptch.xe(i).globalId,ptch.xc(i).globalId) = w{i}.f'.*E + Eglobal(ptch.xe(i).globalId,ptch.xc(i).globalId);
            for d = 1:dim
                Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId) =  2.*w{i}.grad{d}'.*B{d} + Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId);
            end
            Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId) = w{i}.L'.*E + w{i}.f'.*L + Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId);
        end
        L = Lglobal(dataY.inner,:);
        B = Eglobal(dataY.bnd,:);
    end
end

