classdef LQNS < SolverLQNS
    % LQNS - Alias for SolverLQNS (Layered Queueing Network Solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = LQNS(model, varargin)
            % LQNS(MODEL, VARARGIN)
            self@SolverLQNS(model, varargin{:});
        end
    end
end
