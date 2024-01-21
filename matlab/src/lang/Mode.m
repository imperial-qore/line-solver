classdef Mode < NetworkElement
    % An abstract class for a firing mode
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.

    properties
        index;
    end

    methods (Hidden)
        %Constructor
        function self = Mode(name, id)
            % SELF = MODE(NAME)

            self@NetworkElement(name);
            self.index = id;
        end
    end

    methods (Access=public)
        function ind = subsindex(self)
            % IND = SUBSINDEX()

            ind = double(self.index)-1; % 0 based
        end
    end

end
