classdef EventType < Copyable
    % Types of events 
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    % event major classification
    properties (Constant)
        INIT = -1; % model is initialized (time t=0)
        LOCAL = 0; % dummy event
        ARV = 1; % job arrival
        DEP = 2; % job departure
        PHASE = 3; % service advances to next phase, without departure
        READ = 4; % read cache item
        STAGE = 5; % random environment stage change
        ENABLE = 6; % enable mode
        FIRE = 7; % fire mode
        PRE = 8; % consume from a place or queue buffer (no side-effects on server)
        POST = 9; % produce to a place or queue buffer
        RENEGE = 10; % a waiting job abandons the queue (impatience)
        RETRY = 11; % an orbiting job retries entry into a retrial station
        SWITCH = 12; % the server of a polling station advances its switchover timer
        FAILURE = 13; % the server of a station breaks down (goes from up to down)
        REPAIR = 14; % the server of a station is repaired (goes from down to up)
        START = 15; % a job begins or resumes holding a server (tag on an ARV/DEP arc)
        PREEMPT = 16; % a job holding a server is pushed back into the buffer (tag on an ARV arc)
    end
    % START and PREEMPT are instantaneous tags on the arc of the transition
    % that causes them, never the active half of an sn.sync entry: they carry
    % no clock, add no state and leave every numerical result unchanged.
    % PREEMPT is spelled in full because PRE = 8 already names the Petri-net
    % pre-arc. REPAIR emits no START on purpose: the supported breakdown model
    % resumes the held job (downServiceRates degrades, it does not evict).
    
    methods(Static)
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            
            switch type
                case EventType.ARV
                    text = 'ARV';
                case EventType.DEP
                    text = 'DEP';
                case EventType.PHASE
                    text = 'PHASE';
                case EventType.READ
                    text = 'READ';
                case EventType.LOCAL
                    text = 'LOCAL';
                case EventType.STAGE
                    text = 'STAGE';
                case EventType.ENABLE
                    text = 'ENABLE';
                case EventType.FIRE
                    text = 'FIRE';
                case EventType.PRE
                    text = 'PRE';
                case EventType.POST
                    text = 'POST';
                case EventType.RENEGE
                    text = 'RENEGE';
                case EventType.RETRY
                    text = 'RETRY';
                case EventType.SWITCH
                    text = 'SWITCH';
                case EventType.FAILURE
                    text = 'FAILURE';
                case EventType.REPAIR
                    text = 'REPAIR';
                case EventType.START
                    text = 'START';
                case EventType.PREEMPT
                    text = 'PREEMPT';
            end
        end        
    end
    
end
