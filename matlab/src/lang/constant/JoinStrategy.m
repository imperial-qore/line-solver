classdef (Sealed) JoinStrategy
    % Enumeration of join strategy types.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        STD = 'Standard Join';
        PARTIAL = 'Partial Join';
        %QUORUm = 'Partial Join';
        %GUARD = 'Partial Join';
    end
    
    methods (Access = private)
        %private so that you can't instatiate.
        function out = JoinStrategy
            % OUT = JOINSTRATEGY
            
        end
    end
    
    methods (Static)
        
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            
            switch type
                case JoinStrategy.STD
                    text = 'Stardard Join';
                case JoinStrategy.PARTIAL
                    text = 'Partial Join';
            end
        end
    end
    
    
end

