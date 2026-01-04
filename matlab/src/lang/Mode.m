classdef Mode < NetworkElement
    % An abstract class for a firing mode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        transition;
        index;
    end

    methods (Hidden)
        %Constructor
        function self = Mode(T, name)
            % SELF = MODE(NAME)

            self@NetworkElement(name);
            self.index = length(T.getModes()) + 1;
            self.transition = T;
        end
    end

    methods (Access=public)
        function ind = subsindex(self)
            % IND = SUBSINDEX()

            ind = double(self.index)-1; % 0 based
        end
    end

end
