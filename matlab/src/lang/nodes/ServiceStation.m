classdef ServiceStation < Station
    % An abstract class for stations with service
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        schedPolicy;
        schedStrategy;
        schedStrategyPar;
        serviceProcess;
    end

    methods (Hidden)
        %Constructor
        function self = ServiceStation(name)
            % SELF = STATION(NAME)

            self@Station(name);
        end
    end

    methods 
        function distrib = getServiceProcess(self, oclass)
            distrib = self.getService{oclass};
        end
    end
end
