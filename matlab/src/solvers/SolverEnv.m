classdef SolverEnv < SolverENV
    % SolverEnv - Deprecated alias for SolverENV (Ensemble environment solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverEnv(model, varargin)
            % SolverEnv(MODEL, VARARGIN)
            self@SolverENV(model, varargin{:});
            warning('LINE:Deprecated', 'The SolverEnv class is deprecated. Please use SolverENV instead.');
        end
    end
end
