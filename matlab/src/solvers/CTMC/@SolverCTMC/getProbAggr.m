function varargout = getProbAggr(self,varargin)
% PNIR = GETPROBAGGR(IST)
%
% Probability of a SPECIFIC per-class job distribution at a station.
% Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for current state.
%
% Compare with getProbMarg: returns total queue-length distribution,
% i.e., P(n total jobs) summed over all class combinations.
%
% Input:
%   ist - Station index or node object
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

function Pnir = getProbAggr_impl(self, ist)
% GETPROBAGGR_IMPL Implementation of GETPROBAGGR; see the wrapper above.


% lang='cpp' answers this from -a prob; a state prior over several rows is
% refused there by name. See CPPLINE.assertSingleState.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    CPPLINE.assertSingleState(self.name, 'getProbAggr', self.model);
    if ~isnumeric(ist), istc = ist.index; else, istc = ist; end
    Pnir = CPPLINE.probEntry(CPPLINE.probAggr(self.name, self.model, self.options), ...
        'ProbAggr', istc, self.name, 'getProbAggr');
    return
end

self.assertPhaseTypeStates('getProbAggr');

if GlobalConstants.DummyMode
    Pnir = NaN;
    return
end

if ~isnumeric(ist) % station object
    ist = ist.index; 
end

sn = self.getStruct;
if nargin<2 %~exist('ist','var')
    line_error(mfilename,'getProb requires to pass a parameter the station of interest.');
end
if ist > sn.nstations
    line_error(mfilename,'Station number exceeds the number of stations in the model.');
end
if ~isfield(self.options,'keep')
    self.options.keep = false;
end
T0 = tic;
sn.state = sn.state;

if isempty(self.result) || ~isfield(self.result,'Prob') || ~isfield(self.result.Prob,'marginal')
    Pnir = solver_ctmc_margaggr(sn, self.options);
    self.result.('solver') = getName(self);
    self.result.Prob.marginal = Pnir;
else
    Pnir = self.result.Prob.marginal;
end
runtime = toc(T0);
self.result.runtime = runtime;
Pnir = Pnir(ist);
end