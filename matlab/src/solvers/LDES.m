classdef LDES < SolverLDES
    % LDES - Alias for SolverLDES
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = LDES(model, varargin)
            % LDES(MODEL, VARARGIN)
            self@SolverLDES(model, varargin{:});
        end
    end
end
