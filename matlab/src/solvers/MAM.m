classdef MAM < SolverMAM
    % MAM - Alias for SolverMAM (Matrix Analytic Methods solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = MAM(model, varargin)
            % MAM(MODEL, VARARGIN)
            self@SolverMAM(model, varargin{:});
        end
    end
end
