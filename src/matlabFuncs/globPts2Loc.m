%
% Here, we can either have an index list or a structure for xglob
%
function xloc = globPts2Loc(C,T,xglob,plist)
  for k=1:length(T)
    if (nargin==4)  
      if (isstruct(plist))
        locpts = plist(k).ind;
      else
        locpts = find(plist(:,k));
      end  
      xg = xglob(locpts,:); 
    else % Already a structure
      xg = xglob(k).pts;
    end
    xg = xg - ones(size(xg,1),1)*C(k,:);
    xg = xg*T(k).Q;
    
    xloc(k).pts = xg;
  end
  
    
    
    