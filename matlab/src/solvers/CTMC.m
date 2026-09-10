classdef CTMC < SolverCTMC
    % CTMC - Alias for SolverCTMC (Continuous Time Markov Chain solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = CTMC(model, varargin)
            % CTMC(MODEL, VARARGIN)
            self@SolverCTMC(model, varargin{:});
        end
    end
end
