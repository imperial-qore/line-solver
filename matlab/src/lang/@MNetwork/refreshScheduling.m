function [sched, schedparam] = refreshScheduling(self)
% [SCHED, SCHEDPARAM] = REFRESHSCHEDULING()
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% determine scheduling parameters
M = getNumberOfStations(self);
K = getNumberOfClasses(self);

sched = getStationScheduling(self);
schedparam = zeros(M,K);
for ist=1:M
    if isempty(self.getIndexSourceStation) || ist ~= self.getIndexSourceStation
        switch self.stations{ist}.server.className
            case 'ServiceTunnel'
                % do nothing
            otherwise
                if ~isempty(self.stations{ist}.schedStrategyPar) & ~isnan(self.stations{ist}.schedStrategyPar) %#ok<AND2>
                    schedparam(ist,:) = self.stations{ist}.schedStrategyPar;
                else
                    switch sched(ist)
                        case SchedStrategy.SEPT
                            servTime = zeros(1,K);
                            for k=1:K
                                servTime(k) = self.nodes{ist}.serviceProcess{k}.getMean;
                            end
                            [servTimeSorted] = sort(unique(servTime));
                            self.nodes{ist}.schedStrategyPar = zeros(1,K);
                            for k=1:K
                                self.nodes{ist}.schedStrategyPar(k) = find(servTimeSorted == servTime(k));
                            end
                        case SchedStrategy.LEPT
                            servTime = zeros(1,K);
                            for k=1:K
                                servTime(k) = self.nodes{ist}.serviceProcess{k}.getMean;
                            end
                            [servTimeSorted] = sort(unique(servTime),'descend');
                            self.nodes{ist}.schedStrategyPar = zeros(1,K);
                            for k=1:K
                                self.nodes{ist}.schedStrategyPar(k) = find(servTimeSorted == servTime(k));
                            end
                    end
                end
        end
    end
end

if ~isempty(self.sn)
    self.sn.schedparam = schedparam;
	self.sn.sched = sched(:);
end
end
