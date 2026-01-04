function sampleResult = sampleAggr(self, node, numSamples)
% SAMPLERESULT = SAMPLEAGGR(NODE, NUMSAMPLES) Returns aggregated sample path at a node
%
% @brief Generates a discrete-event sample path with state aggregated over phases
%
% This method runs a DES transient simulation and records the sequence
% of states visited at a specified station, with states aggregated to
% show only the number of jobs per class (not per-phase details).
%
% @param self SolverDES instance
% @param node Queue or node object to sample
% @param numSamples Number of time points to sample
%
% @return sampleResult Structure containing the aggregated sample path with fields:
%         - sampleResult.handle: Reference to the sampled node
%         - sampleResult.t: Vector of time points
%         - sampleResult.state: Matrix of aggregated states (jobs per class)
%         - sampleResult.event: Cell array of events
%         - sampleResult.isaggregate: Boolean true (aggregated sample)
%
% @see sample - Returns detailed sample path with phase information
% @see sampleSysAggr - Returns aggregated system-wide sample path
%
% Example:
% @code
% model = Network('MM1');
% % ... model setup ...
% solver = SolverDES(model, 'seed', 42);
% queue = model.nodes{1};
% result = solver.sampleAggr(queue, 1000);
% @endcode

if nargin < 2
    line_error(mfilename, 'sampleAggr requires to specify a station.');
end

if nargin < 3
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

% Get stateful node from model
sn = self.model.getStruct;
if isa(node, 'Node')
    nodeIdx = node.index;
else
    nodeIdx = node;
end

% Find the corresponding Java node
statefulNodes = jmodel.getStatefulNodes();
isf = sn.nodeToStateful(nodeIdx);
jnode = statefulNodes.get(isf - 1);  % Java is 0-indexed

% Call Java sampleAggr method
jresult = jsolver.sampleAggr(jnode, numSamples);

% Convert Java result to MATLAB structure
sampleResult = struct();
sampleResult.handle = node;

% Convert time matrix
if ~isempty(jresult.t)
    sampleResult.t = JLINE.from_jline_matrix(jresult.t);
else
    sampleResult.t = [];
end

% Convert state matrix
if ~isempty(jresult.state) && isa(jresult.state, 'jline.util.matrix.Matrix')
    sampleResult.state = JLINE.from_jline_matrix(jresult.state);
else
    sampleResult.state = [];
end

% Convert event matrix
if ~isempty(jresult.event) && isa(jresult.event, 'jline.util.matrix.Matrix')
    sampleResult.event = JLINE.from_jline_matrix(jresult.event);
else
    sampleResult.event = {};
end

sampleResult.isaggregate = jresult.isAggregate;
sampleResult.nodeIndex = nodeIdx;
sampleResult.numSamples = jresult.numSamples;

end
