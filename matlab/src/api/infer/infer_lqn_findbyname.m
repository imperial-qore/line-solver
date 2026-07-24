function el = infer_lqn_findbyname(container, name)
% INFER_LQN_FINDBYNAME Locate an LQN element by name in a cell container.
%
%   EL = INFER_LQN_FINDBYNAME(CONTAINER, NAME) returns the first element of
%   the cell array CONTAINER whose .name equals NAME, or [] if none matches.
%   Used by the LQN parameter identification helpers to resolve activities and
%   tasks referenced by name in a parameter specification (see INFER_LQN).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

el = [];
for k = 1:numel(container)
    if strcmp(container{k}.name, name)
        el = container{k};
        return
    end
end
end
