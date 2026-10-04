classdef (Sealed) JoinStrategy
    % Enumeration of join strategy types.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        STD = 1;
        PARTIAL = 2;
    end
    
    methods (Static)
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            switch type
                case JoinStrategy.STD
                    text = 'Standard Join';
                case JoinStrategy.PARTIAL
                    text = 'Partial Join';
                otherwise
                    line_error(mfilename, 'Unrecognized join strategy type.');
            end
        end
    end

end
