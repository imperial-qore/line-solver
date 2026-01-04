function sampleAggr = sampleAggr(self, node, numSamples)
% S = SAMPLEAGGR(NODE, NUMSAMPLES) Returns sample path aggregated over service phases
%
% @brief Generates a discrete-event sample path showing job counts per class at a node
%
% This method produces a sample path similar to sample(), but with states aggregated
% over the phases of service distributions. The result shows only the number of jobs
% of each class at the station, without distinguishing internal phase information.
%
% Aggregation is useful for higher-level system analysis when phase information is
% not needed, providing a cleaner view of job class populations at each station.
%
% @param self SolverCTMC instance
% @param node Queue or node object (or node index as integer) to sample
% @param numSamples Integer number of events/transitions to simulate (e.g., 1000)
%
% @return sampleAggr Structure containing the aggregated sample path with fields:
%         - sampleAggr.handle: Reference to the sampled node
%         - sampleAggr.t: Vector of time points when state changes occurred
%         - sampleAggr.state: Matrix of aggregated states at each time point
%           Shape: [numSamples x num_classes] (no phase information)
%         - sampleAggr.event: Cell array of events that triggered transitions
%         - sampleAggr.isaggregate: Boolean flag (true, indicating aggregated)
%         - sampleAggr.nodeIndex: Index of the sampled node
%         - sampleAggr.numSamples: Number of samples collected
%
% @note Aggregated states contain only job counts per class. For detailed phase
%       information, use sample() instead. Supported by CTMC, JMT, and SSA solvers.
%
% @warning Results depend on random seed. Use seed option in constructor for reproducibility.
%
% @see sample - Returns phase-detailed sample path
% @see sample_sys_aggr - Returns aggregated system-wide sample path
%
% Example:
% @code
% model = Network('MM1');
% % ... model setup ...
% solver = SolverCTMC(model, 'seed', 42);
% queue = model.nodes{1};
%
% % Get aggregated sample path with 500 transitions
% result = solver.sampleAggr(queue, 500);
%
% % Access aggregated state sequence
% times = result.t;
% job_counts = result.state;  % [n_class0, n_class1, ...]
%
% % Plot results
% plot(times, job_counts(:,1), 'b-');  % Class 0 job count over time
% hold on;
% plot(times, job_counts(:,2), 'r-');  % Class 1 job count over time
% @endcode

options = self.getOptions;
options.force = true;
if isempty(self.result) || ~isfield(self.result,'infGen')
    runAnalyzer(self);
end
[infGen, eventFilt] = getGenerator(self);
stateSpace = getStateSpace(self);

initState = sn.state;
nst = cumsum([1,cellfun(@length,initState)']);
s0 = cell2mat(initState(:)');

% set initial state
pi0 = zeros(1,size(stateSpace,1));
pi0(matchrow(stateSpace,s0))=1;

% filter all CTMC events as a marked Markovian arrival process
D1 = cellsum(eventFilt);
D0 = infGen-D1;
MMAP = mmap_normalize([{D0},{D1},eventFilt(:)']);

% now sampel the MMAP
[sjt,event,~,~,sts] = mmap_sample(MMAP,numSamples, pi0);

sn = self.getStruct;
sampleAggr = struct();
sampleAggr.handle = node;
sampleAggr.t = cumsum([0,sjt(1:end-1)']');
ind = node.index;
isf = sn.nodeToStateful(ind);
sampleAggr.state = stateSpace(sts,(nst(isf):nst(isf+1)-1));
[~,sampleAggr.state] = State.toMarginal(sn,sn.statefulToNode(isf),sampleAggr.state);

sampleAggr.event = {};
%nodeEvent = false(length(event),1);
%nodeTS = zeros(length(event),1);
for e = 1:length(event)
    for a=1:length(sn.sync{event(e)}.active)
        sampleAggr.event{end+1} = sn.sync{event(e)}.active{a};
        sampleAggr.event{end}.t = sampleAggr.t(e);
%        if  sn.sync{event(e)}.active{a}.node == ind 
%            nodeEvent(e) = true;
%            nodeTS(e) = tranSysState.t(e);
%        end
    end
    for p=1:length(sn.sync{event(e)}.passive)
        sampleAggr.event{end+1} = sn.sync{event(e)}.passive{p};
        sampleAggr.event{end}.t = sampleAggr.t(e);
%        if  sn.sync{event(e)}.passive{p}.node == ind 
%            nodeEvent(e) = true;
%            nodeTS(e) = tranSysState.t(e);
%        end
    end
end

%tranSysState.state = tranSysState.state([1;find(nodeTS>0)],:);
%tranSysState.t = unique(nodeTS);
%tranSysState.event = tranSysState.event(nodeEvent)';
sampleAggr.isaggregate = true;

end