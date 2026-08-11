classdef SampledMetric < Copyable
    % Observed data for a metric
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        t;
        type;
        data;
        node;
        class;
        cond;
        format; % 'timeseries' (default) or 'trace' (per-request)
    end

    methods

        function self = SampledMetric(type, ts, data, node, jobclass)
            self.t = ts;
            self.data = data;
            self.type = type;
            self.node = node;
            self.class = [];
            self.cond = [];
            self.format = 'timeseries';
            if exist('jobclass','var')
                self.class = jobclass;
            end
        end

        function setConditional(self, event)
            if exist('event','var')
                self.cond = event;
            end
        end

        function setTrace(self)
            self.format = 'trace';
        end

        function bool = isAggregate(self)
            bool = isempty(self.class);
        end

        function bool = isConditional(self)
            bool = ~isempty(self.cond);
        end

        function bool = isTrace(self)
            bool = strcmp(self.format, 'trace');
        end
    end
end

