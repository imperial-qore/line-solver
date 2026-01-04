function tranSysState = sampleSysAggr(self, numSamples)
% TRANSYSSTATE = SAMPLESYSAGGR(NUMSAMPLES) Returns aggregated system-wide sample path
%
% @brief Generates a discrete-event sample path for all stateful nodes with aggregation
%
% This method runs a DES transient simulation and records the sequence
% of states visited at all stateful nodes, with states aggregated to
% show only the number of jobs per class (not per-phase details).
%
% @param self SolverDES instance
% @param numSamples Number of time points to sample
%
% @return tranSysState Structure containing the aggregated system-wide sample path with fields:
%         - tranSysState.handle: Cell array of stateful node handles
%         - tranSysState.t: Vector of time points
%         - tranSysState.state: Cell array of aggregated state matrices (jobs per class)
%         - tranSysState.event: Cell array of events
%         - tranSysState.isaggregate: Boolean true (aggregated sample)
%
% @see sampleSys - Returns detailed system-wide sample path
% @see sampleAggr - Returns aggregated sample path for a single node
%
% Example:
% @code
% model = Network('MM1');
% % ... model setup ...
% solver = SolverDES(model, 'seed', 42);
% result = solver.sampleSysAggr(1000);
% @endcode

if nargin < 2
    numSamples = self.options.samples;
end

% Convert model to Java
jmodel = LINE2JLINE(self.model);

% Create Java DES options
joptions = jline.solvers.des.DESOptions();
joptions.samples = self.options.samples;
joptions.seed = self.options.seed;

% Create Java solver
jsolver = jline.solvers.des.SolverDES(jmodel, joptions);

% Call Java sampleSysAggr method
jresult = jsolver.sampleSysAggr(numSamples);

% Convert Java result to MATLAB structure
tranSysState = struct();
tranSysState.handle = self.model.getStatefulNodes';

% Convert time matrix
if ~isempty(jresult.t)
    tranSysState.t = JLINE.from_jline_matrix(jresult.t);
else
    tranSysState.t = [];
end

% Convert state - it's a List<Matrix> in Java, convert to cell array
sn = self.model.getStruct;
if ~isempty(jresult.state) && isa(jresult.state, 'java.util.List')
    nstateful = jresult.state.size();
    tranSysState.state = cell(1, nstateful);
    for isf = 1:nstateful
        jmatrix = jresult.state.get(isf - 1);  % Java is 0-indexed
        if ~isempty(jmatrix) && isa(jmatrix, 'jline.util.matrix.Matrix')
            tranSysState.state{isf} = JLINE.from_jline_matrix(jmatrix);
        else
            tranSysState.state{isf} = [];
        end
    end
else
    tranSysState.state = {};
end

% Convert event matrix
if ~isempty(jresult.event) && isa(jresult.event, 'jline.util.matrix.Matrix')
    tranSysState.event = JLINE.from_jline_matrix(jresult.event);
else
    tranSysState.event = {};
end

tranSysState.isaggregate = jresult.isAggregate;
tranSysState.numSamples = jresult.numSamples;

end
