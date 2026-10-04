classdef FunctionTask < SetupTask
    % FunctionTask is the former name of SetupTask, kept for backward compatibility.
    %
    % Setup and delay-off times are not specific to serverless
    % (function-as-a-service) platforms, so the class carrying them is now
    % named after the modelling primitive rather than after that application
    % domain. Use SetupTask, or a plain Task with setSetupTime() and
    % setDelayOffTime().
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        %constructor
        function self = FunctionTask(model, name, multiplicity, scheduling)
            % self = FunctionTask(model, name, multiplicity, scheduling)
            % Deprecated alias of SetupTask.

            if nargin < 2
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end

            if nargin < 3
                multiplicity = 1;
            end

            if nargin < 4
                scheduling = SchedStrategy.FCFS;
            end

            self@SetupTask(model, name, multiplicity, scheduling);
        end
    end
end
