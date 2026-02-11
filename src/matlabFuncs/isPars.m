% -------------------------------------------------------------------------
% isPars.m      -- Check if input is a structure with the same fields and
%                  types as given by setPars
% Inputs        -- inPars  -> Structure to be checked.
% Outputs       -- outFlag -> Boolean flag.
% Syntax        -- pars.dim = 2; pars.mode = "unfitted"; ...
%                  if isPars(pars) disp("Manual pars is of correct type"); end                   
% Copyright (c) 2025 Andreas Michael <andreas.michael@it.uu.se>
%
% All rights reserved. Use of this source code is governed by a
% BSD-style license that can be found in the LICENSE file.
% -------------------------------------------------------------------------
function [outFlag] = isPars(inPars)
    checkPars = setPars;
    checkPars = orderfields(checkPars);
    inPars = orderfields(inPars);
    flds = fieldnames(checkPars)';
    % Do all fields match
    outFlag = all(ismember(fieldnames(inPars),flds));
    if outFlag == 0
        return;
    end
    % Do all field types match
    values = cellfun(@(f) [checkPars.(f)], flds, 'UniformOutput', false);
    for i = 1:length(values)
        if ~isa(inPars.(flds{i}), class(values{i}))
            outFlag = 0;
            return;
        end
    end
end