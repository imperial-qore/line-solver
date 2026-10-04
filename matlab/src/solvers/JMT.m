classdef JMT < SolverJMT
    % JMT - Alias for SolverJMT (Java Modelling Tools solver)
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = JMT(model, varargin)
            % JMT(MODEL, VARARGIN)
            self@SolverJMT(model, varargin{:});
        end
    end
end
