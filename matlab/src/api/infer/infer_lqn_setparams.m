function model = infer_lqn_setparams(model, paramSpec, a)
% INFER_LQN_SETPARAMS Apply a parameter vector to a LayeredNetwork.
%
%   MODEL = INFER_LQN_SETPARAMS(MODEL, PARAMSPEC, A) sets the LQN parameters
%   named in PARAMSPEC to the values in the vector A, in place, and returns
%   the (same handle) MODEL. This is the parameter-injection half of the LQN
%   parameter identification method (see INFER_LQN).
%
%   PARAMSPEC is a struct array; PARAMSPEC(i) has fields:
%       .type : 'hostdem' (activity mean host demand) or
%               'think'   (task mean think time)
%       .name : name of the activity ('hostdem') or task ('think')
%
%   Mean values are injected as exponential distributions (SCV = 1), which is
%   the assumption of the CASCON 2005 tracking method. The cached layered
%   struct (MODEL.lsn) is invalidated so the next getStruct/solve re-reads the
%   mutated model objects.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

np = numel(paramSpec);
if numel(a) ~= np
    line_error(mfilename, 'Length of parameter vector does not match paramSpec.');
end

for i = 1:np
    val = a(i);
    switch lower(paramSpec(i).type)
        case 'hostdem'
            act = infer_lqn_findbyname(model.activities, paramSpec(i).name);
            if isempty(act)
                line_error(mfilename, sprintf('Activity ''%s'' not found.', paramSpec(i).name));
            end
            act.setHostDemand(val);
        case 'think'
            tsk = infer_lqn_findbyname(model.tasks, paramSpec(i).name);
            if isempty(tsk)
                line_error(mfilename, sprintf('Task ''%s'' not found.', paramSpec(i).name));
            end
            tsk.setThinkTime(val);
        otherwise
            line_error(mfilename, sprintf('Unknown parameter type ''%s''.', paramSpec(i).type));
    end
end

% invalidate the cached layered struct so the next solve re-reads the objects
model.lsn = [];
end
