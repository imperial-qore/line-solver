classdef QNS < SolverQNS
    % QNS - Alias for SolverQNS (Queueing Network Solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = QNS(model, varargin)
            % QNS(MODEL, VARARGIN)
            self@SolverQNS(model, varargin{:});
        end
    end
end
