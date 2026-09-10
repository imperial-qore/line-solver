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

        function self = removeJobClass(self, jobclass)
            % SELF = REMOVEJOBCLASS(JOBCLASS)
            %
            % Drop the service process and the scheduling parameter of
            % JOBCLASS, at the station and inside its server section, on top
            % of what Station and Node remove.

            remaining = self.remainingClassIndexes(jobclass);
            removeJobClass@Station(self, jobclass);
            K = numel(remaining) + 1;
            if numel(self.schedStrategyPar) == K
                self.schedStrategyPar = self.schedStrategyPar(remaining);
            end
            if numel(self.serviceProcess) == K
                self.serviceProcess = self.serviceProcess(remaining);
            end
            if ~isempty(self.server) && isprop(self.server, 'serviceProcess') ...
                    && numel(self.server.serviceProcess) == K
                self.server.serviceProcess = self.server.serviceProcess(remaining);
            end
        end
    end
end
