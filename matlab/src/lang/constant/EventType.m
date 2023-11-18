classdef EventType < Copyable
    % Types of events 
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    % event major classification
    properties (Constant)
        ID_INIT = -1; % model is initialized (t=0)
        ID_LOCAL = 0;
        ID_ARV = 1; % job arrival
        ID_DEP = 2; % job departure
        ID_PHASE = 3; % service advances to next phase, without departure
        ID_READ = 4; % read cache item
        ID_STAGE = 5; % random environment stage change
        
        INIT = 'init'; 
        LOCAL = 'local';
        ARV = 'arv';
        DEP = 'dep';
        PHASE = 'phase';
        READ = 'read';
        STAGE = 'stage';
    end
    
    methods(Static)
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            
            switch type
                case EventType.ID_ARV
                    text = EventType.ARV;
                case EventType.ID_DEP
                    text = EventType.DEP;
                case EventType.ID_PHASE
                    text = EventType.PHASE;
                case EventType.ID_READ
                    text = EventType.READ;
                case EventType.ID_STAGE
                    text = EventType.STAGE;
            end
        end
        function id = toId(type)
            % id = TOID(TYPE)
            
            switch type
                case EventType.ARV
                    id = EventType.ID_ARV;
                case EventType.DEP
                    id = EventType.ID_DEP;
                case EventType.PHASE
                    id = EventType.ID_PHASE;
                case EventType.READ
                    id = EventType.ID_READ;
                case EventType.STAGE
                    id = EventType.ID_STAGE;
            end
        end
    end
    
end
