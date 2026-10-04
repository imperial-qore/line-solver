function varargout = getProb(self,varargin)
% PNIR = GETPROB(NODE, STATE)
% The result recorder captures the scalar this getter returned together
% with the solver that produced it -- see LineResultRecorder. Six of the
% statepr_* goldens hold exactly this number and nothing else, so recording
% it is what makes those goldens attributable instead of "the first bare
% number the example printed". The wrapper exists so that recording happens
% on EVERY exit path of the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getProb_impl(self,varargin{:});
LineResultRecorder.captureScalar(scope, self, 'prob', varargout{1});
end

function Pnir = getProb_impl(self, node, state)
% GETPROB_IMPL Implementation of GETPROB; see the wrapper above.


if GlobalConstants.DummyMode
    Pnir = NaN;
    return
end

if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    % The state, when one is given, crosses as --node/--state and is decoded
    % there with the port's to_marginal, which is the reference's own
    % `sn.state{ist} = state` followed by State.toMarginal. Without one the
    % model's declared state stands, and that already travels in model.json.
    sn = self.getStruct;
    if isa(node,'Node')
        ind = node.index;
    else
        ind = sn.stationToNode(node);
    end
    ist = sn.nodeToStation(ind);
    given = [];
    if nargin >= 3
        given = state;
    end
    p = CPPLINE.probAggr(self.name, self.model, self.options, CPPLINE.stateFlags(ind, given));
    % THIS GETTER RETURNS A LOG, and the name does not say so. Both references
    % agree on it: solver_nc_marg's first output is `lPr` and the native path
    % below returns it unexponentiated, as does SolverNC.java (`Pnir =
    % ret.lPr; return Pnir.get(ist)`). The wire keeps `Prob` a probability,
    % because that key means the same thing under every -s, so the logarithm is
    % taken HERE rather than the CLI reporting a log under a name that says
    % probability everywhere else.
    Pnir = log(CPPLINE.probEntry(p, 'Prob', ist, self.name, 'getProb'));
    return
end

T0 = tic;
% sn was READ THREE LINES ABOVE its assignment, so getProb(node) without an
% explicit state threw "Unrecognized function or variable 'sn'" and only the
% two-argument call ever ran. Same defect the sampleAggr accessors carried.
sn = self.getStruct;
if nargin<3 %~exist('state','var')
    state = sn.state{sn.nodeToStateful(node.index)};
end
% now compute marginal probability
if isa(node,'Node')
    ist = sn.nodeToStation(node.index);
else
    ist = node;    
end
sn.state{ist} = state;

options = self.getOptions;
Solver.resetRandomGeneratorSeed(options.seed);

if ~isempty(self.result) && ~isempty(self.result.Prob) && isfield(self.result.Prob,'logNormConstAggr') && isfinite(self.result.Prob.logNormConstAggr)
    [Pnir,lG] = solver_nc_marg(sn, self.options, self.result.Prob.logNormConstAggr);
else
    [Pnir,lG] = solver_nc_marg(sn, self.options);
    self.result.Prob.logNormConstAggr = lG;
end
self.result.('solver') = getName(self);
self.result.Prob.marginal = Pnir;
runtime = toc(T0);
self.result.runtime = runtime;
Pnir = Pnir(ist);
end