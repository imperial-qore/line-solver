classdef AG < SolverAG
    % AG - Alias for SolverAG (agent-based RCAT solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = AG(model, varargin)
            % AG(MODEL, VARARGIN)
            self@SolverAG(model, varargin{:});
        end
    end
end
