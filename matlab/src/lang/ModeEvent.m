classdef ModeEvent
    % A mode event occurring in a Network.
    %
    % Object of the Event class are not passed by handle.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        node;
        event;
        mode;
        weight;
        prob;
        state; % state information when the event occurs (optional)
        t; % timestamp when the event occurs (optional)
        job; % job id (optional) 
    end
    
    methods
        function self = ModeEvent(event, node, mode, weight, prob, state, t, job)
            % SELF = MODEEVENT(EVENT, NODE, MODE, WEIGHT, PROB, STATE, TIMESTAMP, JOB)
                                    
            self.node = node;            
            self.event = event;
            self.mode = mode;
            if nargin <4
                weight = 1;
            end
            self.weight = weight;
            if nargin <5
                prob = NaN;
            end
            self.prob = prob;
            if nargin <6
                state = []; % local state of the node or environment transition
            end
            self.state = state;
            if nargin <7
                t = NaN; % timestamp
            end
            self.t = t;
            if nargin <8
                job = NaN; % timestamp
            end
            self.job = job;
        end
       
        function print(self)
            % PRINT()
            if isnan(self.t)
                line_printf('(%s: node: %d, mode: %d)',EventType.toText(self.event),self.node,self.mode);
            else
                if isnan(self.job)
                    line_printf('(%s: node: %d, mode: %d, time: %d)',EventType.toText(self.event),self.node,self.mode,self.t);
                else
                    line_printf('(%s: node: %d, mode: %d, job: %d, time: %d)',EventType.toText(self.event),self.node,self.mode,self.job,self.t);
                end
            end
        end
    end
    
    
end
