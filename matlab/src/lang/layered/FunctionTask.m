classdef FunctionTask < Task
    % FunctionTask is an alias for Task, provided for backward compatibility.
    %
    % All setup/delayoff functionality has been moved to the base Task class.
    % Any Task can now have setup time (cold start delay) and delay-off time
    % (teardown delay) configured via setSetupTime() and setDelayOffTime().
    %
    % This class is retained for backward compatibility with existing code
    % that uses FunctionTask to model serverless functions or tasks with
    % initialization overhead.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Dependent)
        % Backward-compatible property names (mapped to parent's camelCase properties)
        SetupTime;
        SetupTimeMean;
        SetupTimeSCV;
        DelayOffTime;
        DelayOffTimeMean;
        DelayOffTimeSCV;
    end

    methods
        %constructor
        function self = FunctionTask(model, name, multiplicity, scheduling)
            % self = FunctionTask(model, name, multiplicity, scheduling)
            % Create a FunctionTask (alias for Task with setup/delayoff support).

            if nargin < 2
                line_error(mfilename,'Constructor requires to specify at least a name.');
            end

            if nargin < 3
                multiplicity = 1;
            end

            if nargin < 4
                scheduling = SchedStrategy.FCFS;
            end

            self@Task(model, name, multiplicity, scheduling);
        end

        % Dependent property getters/setters for backward compatibility
        function value = get.SetupTime(self)
            value = self.setupTime;
        end

        function set.SetupTime(self, value)
            self.setSetupTime(value);
        end

        function value = get.SetupTimeMean(self)
            value = self.setupTimeMean;
        end

        function set.SetupTimeMean(self, value)
            self.setupTimeMean = value;
        end

        function value = get.SetupTimeSCV(self)
            value = self.setupTimeSCV;
        end

        function set.SetupTimeSCV(self, value)
            self.setupTimeSCV = value;
        end

        function value = get.DelayOffTime(self)
            value = self.delayOffTime;
        end

        function set.DelayOffTime(self, value)
            self.setDelayOffTime(value);
        end

        function value = get.DelayOffTimeMean(self)
            value = self.delayOffTimeMean;
        end

        function set.DelayOffTimeMean(self, value)
            self.delayOffTimeMean = value;
        end

        function value = get.DelayOffTimeSCV(self)
            value = self.delayOffTimeSCV;
        end

        function set.DelayOffTimeSCV(self, value)
            self.delayOffTimeSCV = value;
        end
    end
end
