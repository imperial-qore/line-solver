function sched = getStationScheduling(self)
% SCHED = GETSTATIONSCHEDULING()

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.
for i=1:getNumberOfStations(self)
    if isinf(self.stations{i}.numberOfServers)
        sched{i,1} = SchedStrategy.INF;
    else
        if i == getIndexSourceStation(self)
            sched{i,1} = SchedStrategy.EXT;
        else
            sched{i,1} = self.stations{i}.schedStrategy;
        end
    end
end
end
