classdef (Sealed) ParameterType
    % An input parameter of a Solver
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        ServDemand = 'Service Demand'; % Service Time * Visits
        ServTime = 'Service Time'; % Service Time * Visits
        Visits = 'Visits'; % Service Time * Visits
        NumJobs = 'Number of Jobs'; % Closed Population
        ArvRate = 'Arrival Rate'; % Class Arrival Rate
        ThinkTime = 'Think Time'; % Class Arrival Rate
        
        ID_ServDemand = 0;
        ID_ServTime = 1;
        ID_Visits = 2;
        ID_NumJobs = 3;
        ID_ArvRate = 4;
        ID_ThinkTime = 5;
    end
        
end

