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
    end
    
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
            end
        end        
    end
    
end
