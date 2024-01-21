classdef closedWorkload
    % Auxiliary data structure for closed classes in PMIF models
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties
        name;               %string
        numberJobs;         %int
        thinkTime;          %string
        thinkDevice;        %double
        timeUnits = '';     %string - optional
        transits;           %cell
    end
    
    methods
        %public methods, including constructor
        
        %constructor
        function obj = closedWorkload(name, numberJobs, thinkTime, thinkDevice, timeUnits)
            % OBJ = CLOSEDWORKLOAD(NAME, NUMBERJOBS, THINKTIME, THINKDEVICE, TIMEUNITS)
            
            if(nargin > 0)
                obj.name = name;
                obj.numberJobs = numberJobs;
                obj.thinkTime = thinkTime;
                obj.thinkDevice = thinkDevice;
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
