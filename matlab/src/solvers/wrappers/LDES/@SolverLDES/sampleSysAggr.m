function tranSysState = sampleSysAggr(self, numEvents)
% TRANSYSSTATE = SAMPLESYSAGGR(NUMEVENTS)
% Aggregated system-wide sample path. For LDES the trajectories are already
% per-class, so this equals sampleSys() marked aggregated. Fully JSON-mediated.
% Mirrors the Python-native sampleSysAggr().

if nargin < 2 || isempty(numEvents)
    numEvents = 0;
end
tranSysState = self.sampleSys(numEvents);
if ~isempty(tranSysState) && isstruct(tranSysState)
    tranSysState.isaggregate = true;
end
end
