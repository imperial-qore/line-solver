classdef SetupTask < Task
    % SetupTask is a Task whose servers are switched off while idle.
    %
    % A server resuming from the off state pays a setup (activation) time
    % before serving the request that woke it up, and stays available for a
    % delay-off (idle) period after emptying its queue before switching off.
    % These are the setup and close-down times of a server with vacations:
    % on-demand virtual machines and containers, power-managed servers under a
    % timeout policy, warm-up delays, serverless cold start / keep-alive.
    %
    % Both times are declared on the base Task class via setSetupTime() and
    % setDelayOffTime(), so this subclass is a naming convenience.
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
        function self = SetupTask(model, name, multiplicity, scheduling)
            % self = SetupTask(model, name, multiplicity, scheduling)
            % Create a Task that pays setup and delay-off times.

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
