function sched = getStationScheduling(self)
% SCHED = GETSTATIONSCHEDULING()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
for i=1:getNumberOfStations(self)
    if isinf(self.stations{i}.numberOfServers)
        sched(i) = SchedStrategy.INF;
    else
        if i == getIndexSourceStation(self)
            sched(i) = SchedStrategy.EXT;
        else
            sched(i) = self.stations{i}.schedStrategy;
        end
    end
end
end
