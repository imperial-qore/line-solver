function varargout = getProbAggr(self,varargin)
% PNIR = GETPROBAGGR(NODE, STATE_A)
%
% Probability of a SPECIFIC per-class job distribution at a station.
% Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) at the node.
%
% Compare with getProbMarg: returns total queue-length distribution,
% i.e., P(n total jobs) summed over all class combinations.
%
% Input:
%   node    - Queue or node object
%   state_a - Per-class job counts, e.g., [2,1] = 2 class-1, 1 class-2
%
% Output:
%   Pnir - Scalar probability in [0,1]
% The result recorder captures the scalar this getter returned together
% with the solver that produced it -- see LineResultRecorder. Six of the
% statepr_* goldens hold exactly this number and nothing else, so recording
% it is what makes those goldens attributable instead of "the first bare
% number the example printed". The wrapper exists so that recording happens
% on EVERY exit path of the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getProbAggr_impl(self,varargin{:});
LineResultRecorder.captureScalar(scope, self, 'probAggr', varargout{1});
end

function Pnir = getProbAggr_impl(self, node, state_a)
% GETPROBAGGR_IMPL Implementation of GETPROBAGGR; see the wrapper above.


% lang='cpp' answers this from -a prob; a state prior over several rows is
% refused there by name. See CPPLINE.assertSingleState.
%
% ONLY FOR THE MODEL'S OWN STATE. `-a prob` reports every station's marginal at
% the state model.json carries; a state_a named in the CALL would have to be
% written onto the node first, and doing that here would edit the caller's model
% to ask a question about it. So an explicit state_a stays with the native path,
% where the reference substitutes it into a COPY of sn.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    if nargin >= 3
        CPPLINE.cppUnsupported(self.name, 'getProbAggr(node, state_a)', ...
            ['line-cli reports every station''s marginal at the state the model carries and ' ...
            'takes no per-call state, so answering a named state_a would mean writing it ' ...
            'onto the caller''s model first']);
    end
    CPPLINE.assertSingleState(self.name, 'getProbAggr', self.model);
    sncpp = self.model.getStruct();
    Pnir = CPPLINE.probEntry(CPPLINE.probAggr(self.name, self.model, self.options), ...
        'ProbAggr', sncpp.nodeToStation(node.index), self.name, 'getProbAggr');
    return
end

if GlobalConstants.DummyMode
    Pnir = NaN;
    return
end

T0 = tic;
sn = self.model.getStruct(true); % sync node states into sn.state

% Get unnormalized probability using original NC method
ist = sn.nodeToStation(node.index);
isf = sn.nodeToStateful(node.index);
if nargin<3
    state_a = sn.state{isf};
else
    % state_a is a per-class marginal count vector; convert it to the
    % internal state encoding expected by State.toMarginal, and store it
    % under the stateful index used by solver_nc_margaggr
    state_a = State.fromMarginal(sn, node.index, state_a);
end

% Store original state and set requested state
original_state = sn.state{isf};
sn.state{isf} = state_a;

options = self.getOptions;
Solver.resetRandomGeneratorSeed(options.seed);

self.result.('solver') = getName(self);
% note: solver_nc_margaggr returns [Pr, G, lG, runtime]; the log-scale
% normalizing constant is the THIRD output, not the second
if isfield(self.result,'Prob') && isfield(self.result.Prob,'logNormConstAggr') && isfinite(self.result.Prob.logNormConstAggr)
    [Pnir_vec,~,lG] = solver_nc_margaggr(sn, self.options, self.result.Prob.logNormConstAggr);
else
    [Pnir_vec,~,lG] = solver_nc_margaggr(sn, self.options);
    self.result.Prob.logNormConstAggr = lG;
end
self.result.Prob.marginal = Pnir_vec;

% solver_nc_margaggr already returns normalized probabilities
% (it computes P = F_i * G(-i) / G which is properly normalized)
% So we simply extract the probability for this station
Pnir = Pnir_vec(ist);

% Restore original state
sn.state{isf} = original_state;

runtime = toc(T0);
self.result.runtime = runtime;

end
