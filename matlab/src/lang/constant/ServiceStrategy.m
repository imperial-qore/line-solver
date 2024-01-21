classdef (Sealed) ServiceStrategy
    % Enumeration of service strategies
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        LI = 'LoadIndependent';
        LD = 'LoadDependent';
        CD = 'ClassDependent';
        SD = 'StateDependent';
        ID_LI = 1; %LoadIndependent
        ID_LD = 2; %LoadDependent
        ID_CD = 3; %ClassDependent
        ID_SD = 4; %StateDependent
    end
    
    methods (Access = private)
        %private so that you can't instatiate.
        function out = ServiceStrategy
            % OUT = SERVICESTRATEGY
            
        end
    end
    
end

