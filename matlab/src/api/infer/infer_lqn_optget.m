function v = infer_lqn_optget(s, f, d)
% INFER_LQN_OPTGET Return option field s.(f), or default d when absent/empty.
%
%   V = INFER_LQN_OPTGET(S, F, D) returns S.(F) if S is a struct with a
%   non-empty field F, otherwise the default value D. Small helper shared by
%   the LQN parameter identification routines (see INFER_LQN).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end
