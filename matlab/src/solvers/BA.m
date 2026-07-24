classdef BA < SolverBA
    % BA - Alias for SolverBA (Bound Analysis solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = BA(model, varargin)
            % BA(MODEL, VARARGIN)
            self@SolverBA(model, varargin{:});
        end
    end
end
