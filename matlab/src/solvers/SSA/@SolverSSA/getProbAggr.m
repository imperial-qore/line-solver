function varargout = getProbAggr(self,varargin)
% PROBAGGR = GETPROBAGGR(NODE, STATE)
%
% Probability of a SPECIFIC per-class job distribution at a station.
% Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for given state.
%
% Compare with getProbMarg: returns total queue-length distribution,
% i.e., P(n total jobs) summed over all class combinations.
%
% Input:
%   node  - Node object
%   state - Per-class job counts, e.g., [2,1] = 2 class-1, 1 class-2
%
% Output:
%   ProbAggr - Scalar probability in [0,1] (estimated via simulation)
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

function ProbAggr = getProbAggr_impl(self, node, state)
% GETPROBAGGR_IMPL Implementation of GETPROBAGGR; see the wrapper above.


% lang='cpp' answers this from -a prob, the simulated occupancy of the state the
% model carries; a prior over several rows is refused by name, and a state named
% in the CALL stays with the native path. See CPPLINE.assertSingleState.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    if nargin >= 3
        CPPLINE.cppUnsupported(self.name, 'getProbAggr(node, state)', ...
            ['line-cli reports the occupancy of the state the model carries and takes no ' ...
            'per-call state, so answering a named one would mean writing it onto the ' ...
            'caller''s model first']);
    end
    CPPLINE.assertSingleState(self.name, 'getProbAggr', self.model);
    if ~isnumeric(node), istc = node.index; else, istc = node; end
    snc = self.getStruct;
    istc = snc.nodeToStation(istc);  % the native reads nodeToStateful(node.index)
    ProbAggr = CPPLINE.probEntry(CPPLINE.probAggr(self.name, self.model, self.options), ...
        'ProbAggr', istc, self.name, 'getProbAggr');
    return
end

if GlobalConstants.DummyMode
    ProbAggr = NaN;
    return
end

% we do not use probSysState as that is for joint states
TranSysStateAggr = self.sampleSysAggr;
sn = self.getStruct;
isf = sn.nodeToStateful(node.index);
TSS = cell2mat({TranSysStateAggr.t,TranSysStateAggr.state{isf}});
TSS(:,1)=[TSS(1,1);diff(TSS(:,1))];
if nargin<3 %~exist('state','var')
    state = sn.state{isf};
end
rows = findrows(TSS(:,2:end), state);
if ~isempty(rows)
    ProbAggr = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.\n');
    ProbAggr = 0;
end
end