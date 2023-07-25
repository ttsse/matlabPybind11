% ----------------------------------------------------------------------------------
% p_func.m -- Evaluate symbolic expression f_in, which is discontinuous at points 
%             x_m, at points x.
% Copyright (c) 2023 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% ----------------------------------------------------------------------------------

function [f] = p_func(f_in,x,x_m,cell_sz)

if cell_sz + 1 ~= length(x_m)
    disp("Give break points equal to number of functions - 1")
    return 
end

i_rows = size(f_in{1},1);
j_cols = size(f_in{1},2);

syms x1 x2 x3

for c = length(x_m):-1:2

    intrval = find(x(:,1) >= x_m(c-1) & x(:,1) < x_m(c));

    for i = 1:i_rows
        for j = 1:j_cols
            f{i,j} = matlabFunction(f_in{c-1}(i,j),'Vars',[x1 x2 x3]);
            f_temp(i,j,intrval) = f{i,j}(x(intrval,1),x(intrval,2),x(intrval,3));
        end
    end
end

if i_rows == 1 && j_cols == 1
    f = reshape(f_temp,[size(x,1),1]);
elseif i_rows == 3 && j_cols == 1
    f = [reshape(f_temp(1,1,:),[size(x,1),1]); ...
         reshape(f_temp(2,1,:),[size(x,1),1]); ...
         reshape(f_temp(3,1,:),[size(x,1),1])];
elseif i_rows == 1 && j_cols == 3
    f = [reshape(f_temp(1,1,:),[size(x,1),1]); ...
         reshape(f_temp(1,2,:),[size(x,1),1]); ...
         reshape(f_temp(1,3,:),[size(x,1),1])];    
else 
    f = f_temp;
end

end