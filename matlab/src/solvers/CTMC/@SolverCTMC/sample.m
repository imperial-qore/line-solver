function S = sample(self, node, numSamples)
% S = SAMPLE(NODE, NUMSAMPLES) Returns a sample path of state evolution at a node
%
% @brief Generates a discrete-event sample path showing state transitions at a node
%
% This method produces a sample path (stochastic simulation trace) showing how the
% number of jobs at a specified station evolves over time. The sample path includes
% detailed phase information from the underlying Markov chain, distinguishing between
% different phases of service distributions.
%
% The sample path is obtained by simulating the CTMC underlying the queueing network,
% filtering events to focus on the station of interest, and recording the sequence
% of states visited along with transition times.
%
% @param self SolverCTMC instance
% @param node Queue or node object (or node index as integer) to sample
% @param numSamples Integer number of events/transitions to simulate (e.g., 1000)
%
% @return S Structure containing the sample path with fields:
%         - S.handle: Reference to the sampled node
%         - S.t: Vector of time points when transitions occurred
%         - S.state: Matrix of states at each time point, with phase information
%           Shape: [numSamples x num_classes+1] where last column is phase number
%         - S.event: Cell array of events that triggered transitions
%         - S.isaggregate: Boolean flag (false, indicating phase-detailed)
%         - S.nodeIndex: Index of the sampled node
%         - S.numSamples: Number of samples collected
%
% @note Phase-detailed sample paths include internal state information from
%       phase-type distributions. Use sampleAggr() for simpler class-based states.
%       Only supported by CTMC and SSA solvers.
%
% @warning Results depend on the random seed and CTMC state space.
%          Use seed option in solver constructor for reproducibility.
%
% @see sampleAggr - Returns sample path aggregated over phases
% @see sample_sys - Returns system-wide sample path for all nodes
%
% Example:
% @code
% solver = SolverCTMC(model, 'seed', 12345);
% queue1 = model.nodes{1};
%
% % Get sample path with 1000 transitions
% S = solver.sample(queue1, 1000);
%
% % Plot sample path
% figure;
% plot(S.t, S.state(:,1));  % Plot state for class 1
% xlabel('Time');
% ylabel('Number of Jobs');
% title('Sample Path at Queue1');
% @endcode

options = self.getOptions;
options.force = true;
if isempty(self.result) || ~isfield(self.result,'infGen')
    runAnalyzer(self);
end
[infGen, eventFilt] = getGenerator(self);
stateSpace = getStateSpace(self);
%stateSpaceAggr = getStateSpaceAggr(self);

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
S = struct();
S.handle = node;
S.t = cumsum([0,sjt(1:end-1)']');
ind = node.index;
isf = sn.nodeToStateful(ind);
S.state = stateSpace(sts,(nst(isf):nst(isf+1)-1));

S.event = {};
%nodeEvent = false(length(event),1);
%nodeTS = zeros(length(event),1);
for e = 1:length(event)
    for a=1:length(sn.sync{event(e)}.active)
        S.event{end+1} = sn.sync{event(e)}.active{a};
        S.event{end}.t = S.t(e);
%        if  sn.sync{event(e)}.active{a}.node == ind 
%            nodeEvent(e) = true;
%            nodeTS(e) = tranSysState.t(e);
%        end
    end
    for p=1:length(sn.sync{event(e)}.passive)
        S.event{end+1} = sn.sync{event(e)}.passive{p};
        S.event{end}.t = S.t(e);
%        if  sn.sync{event(e)}.passive{p}.node == ind 
%            nodeEvent(e) = true;
%            nodeTS(e) = tranSysState.t(e);
%        end
    end
end

%tranSysState.state = tranSysState.state([1;find(nodeTS>0)],:);
%tranSysState.t = unique(nodeTS);
%tranSysState.event = tranSysState.event(nodeEvent)';
S.isaggregate = false;

end