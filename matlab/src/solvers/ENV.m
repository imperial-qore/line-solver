classdef ENV < SolverENV
    % ENV - Alias for SolverENV (Ensemble environment solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = ENV(model, varargin)
            % ENV(MODEL, VARARGIN)
            self@SolverENV(model, varargin{:});
        end
    end
end