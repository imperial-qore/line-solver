classdef (Sealed) DropStrategy
    % Enumeration of drop policies in stations.
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        WaitingQueue = 'waitq';
        Drop = 'drop';
        BlockingAfterService = 'bas';
        
        ID_WAITQ = -1;
        ID_DROP = 1;
        ID_BAS = 2;
    end
    
    methods (Static)
        function type = fromId(id)
            switch id
                case DropStrategy.ID_WAITQ
                    type = DropStrategy.WaitingQueue;
                case DropStrategy.ID_DROP
                    type = DropStrategy.Drop;
                case DropStrategy.ID_BAS
                    type = DropStrategy.BlockingAfterService;
            end
        end
        
        function id = toId(type)
            switch type
                case DropStrategy.WaitingQueue
                    id = DropStrategy.ID_WAITQ;
                case DropStrategy.Drop
                    id = DropStrategy.ID_DROP;
                case DropStrategy.BlockingAfterService
                    id = DropStrategy.ID_BAS;
            end
        end
        
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            
            switch type
                case DropStrategy.WaitingQueue
                    text = 'waiting queue';
                case DropStrategy.Drop
                    text = 'drop';
                case DropStrategy.BlockingAfterService
                    text = 'BAS blocking';
            end
        end
    end

end

