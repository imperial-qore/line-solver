classdef demandServiceRequest
    % Auxiliary data structure for service demands in PMIF models
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties
        workloadName;   %string
        serverID;       %string
        serviceDemand;  %double
        numberVisits;   %int
        timeUnits = ''; %string - optional
        transits;
    end
    
    methods
        %public methods, including constructor
        
        %constructor
        function obj = demandServiceRequest(workloadName, serverID, serviceDemand, numberVisits, timeUnits)
            % OBJ = DEMANDSERVICEREQUEST(WORKLOADNAME, SERVERID, SERVICEDEMAND, NUMBERVISITS, TIMEUNITS)
            
            if(nargin > 0)
                obj.workloadName = workloadName;
                obj.serverID = serverID;
                obj.serviceDemand = serviceDemand;
                obj.numberVisits = numberVisits;
                if nargin > 4
                    obj.timeUnits = timeUnits;
                end
            end
        end
        
        function obj = addTransit(obj, dest, prob)
            % OBJ = ADDTRANSIT(OBJ, DEST, PROB)
            
            if isempty(obj.transits)
                obj.transits = cell(1,2);
                obj.transits{1,1} = dest;
                obj.transits{1,2} = prob;
            else
                obj.transits{end+1,1} = dest;
                obj.transits{end,2} = prob;
            end
        end
        
    end
    
end
