function sampleResult = sampleAggr(self, node, numEvents)
% SAMPLERESULT = SAMPLEAGGR(NODE, NUMEVENTS)
% Aggregated sample path for a stateful node. For LDES the sample path is
% already per-class queue lengths, so this equals sample() marked aggregated.
% Fully JSON-mediated (see sample()). Mirrors the Python-native sampleAggr().

if nargin < 3 || isempty(numEvents)
    numEvents = 0;
end
sampleResult = self.sample(node, numEvents);
if ~isempty(sampleResult) && isstruct(sampleResult)
    sampleResult.isaggregate = true;
end
end
