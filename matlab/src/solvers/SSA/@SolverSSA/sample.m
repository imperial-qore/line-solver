function sampleNodeState = sample(self, node, numSamples, markActivePassive)
% TRANNODESTATE = SAMPLE(NODE, NUMSAMPLES, MARKACTIVEPASSIVE)
% Returns a sample path simulation trace showing state evolution at a node
%
% @brief Generates a stochastic simulation sample path with phase information at a node
%
% This method runs a discrete-event stochastic simulation and records the sequence
% of states visited at a specified station along with timing information. The sample
% path includes detailed phase information from the service distributions, capturing
% all internal phases experienced during service.
%
% The SSA solver uses the Next Reaction Method (NRM) for efficient stochastic
% sampling when possible, or general serial CTMC simulation otherwise.
%
% @param self SolverSSA instance
% @param node Queue or node object (or node index as integer) to sample
% @param numSamples (optional) Number of events to simulate (default: solver's samples option)
%                    Use more samples for longer traces and better statistics.
% @param markActivePassive (optional) Boolean flag to distinguish active vs passive events
%                          - false (default): Mixed event list
%                          - true: Separate 'active' and 'passive' event arrays per time
%
% @return sampleNodeState Structure containing the sample path with fields:
%         - sampleNodeState.handle: Reference to the sampled node
%         - sampleNodeState.t: Vector of event times (increasing sequence)
%         - sampleNodeState.state: Matrix of states at event times
%           Shape: [num_events x (num_classes + phase_dimension)]
%         - sampleNodeState.event: Cell array of event information at each time
%         - sampleNodeState.isaggregate: Boolean false (phase-detailed sample)
%         - sampleNodeState.nodeIndex: Index of the sampled node
%         - sampleNodeState.numSamples: Number of events generated
%
% @note Phase information allows detailed analysis of service behavior including
%       time spent in each phase of phase-type distributions. For simpler class-only
%       analysis, use sampleAggr() instead.
%       Supported by CTMC and SSA solvers.
%
% @warning Simulation-based trace - results are stochastic and depend on random seed.
%          Set seed option in solver constructor for reproducibility:
%          solver = SolverSSA(model, 'seed', 12345);
%
% @see sampleAggr - Returns sample path aggregated over phases (simpler output)
% @see sample_sys - Returns system-wide sample path for all stations
%
% Example:
% @code
% model = Network('M_M_1');
% % ... model setup ...
% solver = SolverSSA(model, 'samples', 5000, 'seed', 54321);
% queue1 = model.nodes{1};
%
% % Generate sample path with 1000 events
% S = solver.sample(queue1, 1000);
%
% % Access and plot sample path
% times = S.t;
% states = S.state;
%
% % Plot number of jobs of class 1 vs time
% figure;
% plot(times, states(:,1));
% xlabel('Time');
% ylabel('Number of Jobs (Class 1)');
% title('Sample Path at Queue');
%
% % Mark active/passive events if needed
% S_marked = solver.sample(queue1, 1000, true);
% @endcode

if GlobalConstants.DummyMode
    sampleNodeState = NaN;
    return
end

options = self.getOptions;

if nargin<4
    markActivePassive = false;
end

if nargin>=3 %exist('numSamples','var')
    options.samples = numSamples;
else
    numSamples = options.samples;
end
switch options.method
    case {'default','serial'}
        options.method = 'serial'; % nrm does not support tran*
        [~, tranSystemState, tranSync] = self.runAnalyzer(options);
        event = tranSync;
        sn = self.getStruct;
        isf = sn.nodeToStateful(node.index);
        sampleNodeState = struct();
        sampleNodeState.handle = node;
        sampleNodeState.t = tranSystemState{1};
        sampleNodeState.state = tranSystemState{1+isf};
        
        sn = self.getStruct;
        sampleNodeState.event = {};
        for e = 1:length(event)
            for a=1:length(sn.sync{event(e)}.active)
                sampleNodeState.event{end+1} = sn.sync{event(e)}.active{a};
                sampleNodeState.event{end}.t = sampleNodeState.t(e);
            end
            for p=1:length(sn.sync{event(e)}.passive)
                sampleNodeState.event{end+1} = sn.sync{event(e)}.passive{p};
                sampleNodeState.event{end}.t = sampleNodeState.t(e);
            end
        end
        sampleNodeState.isaggregate = false;

    otherwise
        line_error(mfilename,'sample is not available in SolverSSA with the chosen method.');
end
%sampleNodeState.t = [0; sampleNodeState.t(2:end)];

if markActivePassive
    apevent = cell(1,length(sampleNodeState.t)-1);
    for ti = 1:length(apevent)
        apevent{ti} = struct('active',[],'passive',[]);
    end
    for e=1:length(sampleNodeState.event)
        ti = find(sampleNodeState.event{e}.t == sampleNodeState.t);
        if ~isempty(ti) && ti<length(sampleNodeState.t)
        switch sampleNodeState.event{e}.event
            case EventType.ARV
                apevent{ti}.passive = sampleNodeState.event{e};
            otherwise
                apevent{ti}.active = sampleNodeState.event{e};
        end
        end
    end
    sampleNodeState.event = apevent';
end

end