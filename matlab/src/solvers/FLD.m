classdef FLD < SolverFluid
    % FLD - Alias for SolverFluid (Fluid/Mean-Field Approximation solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = FLD(model, varargin)
            % FLD(MODEL, VARARGIN)
            self@SolverFluid(model, varargin{:});
        end
    end
end
