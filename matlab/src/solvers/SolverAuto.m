classdef SolverAuto < SolverAUTO
    % SolverAuto - Backward compatibility alias for SolverAUTO
    %
    % This class provides backward compatibility with the old lowercase name.
    % New code should use SolverAUTO instead.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverAuto(model, varargin)
            % SolverAuto(MODEL, VARARGIN)
            % Backward compatibility constructor
            self@SolverAUTO(model, varargin{:});
        end
    end
end
