classdef JobClassType < NetworkElement
    % An abstract class for a collection of indistinguishable jobs
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        DISABLED = -1;
        OPEN = 0;
        CLOSED = 1;
    end

    methods (Static)
        function txt = toText(classType)
            switch classType
                case JobClassType.DISABLED
                    txt = 'disabled';
                case JobClassType.OPEN
                    txt = 'open';
                case JobClassType.CLOSED
                    txt = 'closed';
                otherwise
                    line_error(mfilename, 'Unknown job class type.');
            end
        end
    end
end
