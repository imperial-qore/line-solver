classdef NetworkElement < Element
    % A generic element of a Network model.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
     properties (Hidden)
         attribute; % custom attribute
     end
    
    methods
        %Constructor
        function self = NetworkElement(name)
            % SELF = NETWORKELEMENT(NAME)            
            self@Element(name);
        end        
    end
    
end
