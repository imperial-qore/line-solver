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
            % OI (order-independent) is a pass-and-swap specialization with an
            % empty/zero swap graph. Normalize it to PAS in the struct so every
            % solver/state dispatch site (which tests sn.sched == PAS) handles it
            % via the same order-independent machinery.
            if SchedStrategy.toId(sched(i)) == SchedStrategy.OI
                sched(i) = SchedStrategy.PAS;
            end
        end
    end
end
end
