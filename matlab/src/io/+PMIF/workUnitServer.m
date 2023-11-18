classdef workUnitServer
    % Auxiliary data structure for stations in PMIF models
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties
        name;               %string
        quantity = 1;       %int
        scheduling;         %string
        serviceTime;        %double
        timeUnits = '';     %string - optional
    end
    
    methods
        %public methods, including constructor
        
        %constructor
        function obj = workUnitServer(name, quantity, scheduling, serviceTime, timeUnits)
            % OBJ = WORKUNITSERVER(NAME, QUANTITY, SCHEDULING, SERVICETIME, TIMEUNITS)
            
            if(nargin > 0)
                obj.name = name;
                obj.quantity = quantity;
                obj.scheduling = scheduling;
                obj.serviceTime = serviceTime;
                if nargin > 4
                    obj.timeUnits = timeUnits;
                end
            end
        end
        
    end
    
end
