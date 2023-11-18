function [sched, schedid, schedparam] = refreshScheduling(self)
% [SCHED, SCHEDID, SCHEDPARAM] = REFRESHSCHEDULING()
%
% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

% determine scheduling parameters
M = getNumberOfStations(self);
K = getNumberOfClasses(self);

sched = getStationScheduling(self);
schedparam = zeros(M,K);
for i=1:M
    if isempty(self.getIndexSourceStation) || i ~= self.getIndexSourceStation
        switch self.stations{i}.server.className
            case 'ServiceTunnel'
                % do nothing
            otherwise
                if ~isempty(self.stations{i}.schedStrategyPar) & ~isnan(self.stations{i}.schedStrategyPar) %#ok<AND2>
                    schedparam(i,:) = self.stations{i}.schedStrategyPar;
                else
                    switch SchedStrategy.toId(sched{i})
                        case SchedStrategy.ID_SEPT
                            servTime = zeros(1,K);
                            for k=1:K
                                servTime(k) = self.nodes{i}.serviceProcess{k}.getMean;
                            end
                            [servTimeSorted] = sort(unique(servTime));
                            self.nodes{i}.schedStrategyPar = zeros(1,K);
                            for k=1:K
                                self.nodes{i}.schedStrategyPar(k) = find(servTimeSorted == servTime(k));
                            end
                        case SchedStrategy.ID_LEPT
                            servTime = zeros(1,K);
                            for k=1:K
                                servTime(k) = self.nodes{i}.serviceProcess{k}.getMean;
                            end
                            [servTimeSorted] = sort(unique(servTime),'descend');
                            self.nodes{i}.schedStrategyPar = zeros(1,K);
                            for k=1:K
                                self.nodes{i}.schedStrategyPar(k) = find(servTimeSorted == servTime(k));
                            end
                    end
                end
        end
    end
end

if ~isempty(self.sn)
    self.sn.sched = sched;
    self.sn.schedparam = schedparam;
    schedid = zeros(self.sn.nstations,1);
    for i=1:self.sn.nstations
		schedid(i) = SchedStrategy.toId(sched{i});
    end
	self.sn.schedid = schedid;
end
end
