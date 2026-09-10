function varargout = getProbSysAggr(self,varargin)
% PROBSYSSTATEAGGR = GETPROBSYSSTATEAGGR()
% The result recorder captures the scalar this getter returned together
% with the solver that produced it -- see LineResultRecorder. Six of the
% statepr_* goldens hold exactly this number and nothing else, so recording
% it is what makes those goldens attributable instead of "the first bare
% number the example printed". The wrapper exists so that recording happens
% on EVERY exit path of the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getProbSysAggr_impl(self,varargin{:});
LineResultRecorder.captureScalar(scope, self, 'probSysAggr', varargout{1});
end

function ProbSysAggr = getProbSysAggr_impl(self)
% GETPROBSYSAGGR_IMPL Implementation of GETPROBSYSAGGR; see the wrapper above.


% lang='cpp' answers this from -a prob. See CPPLINE.assertSingleState.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    CPPLINE.assertSingleState(self.name, 'getProbSysAggr', self.model);
    ProbSysAggr = CPPLINE.probAggr(self.name, self.model, self.options).ProbSysAggr;
    return
end

if GlobalConstants.DummyMode
    ProbSysAggr = NaN;
    return
end

TranSysStateAggr = self.sampleSysAggr;
TSS = cell2mat([TranSysStateAggr.t,TranSysStateAggr.state(:)']);
TSS(:,1)=[TSS(1,1);diff(TSS(:,1))];
sn = self.getStruct;
state = sn.state;
nir = zeros(sn.nstateful,sn.nclasses);
for isf=1:sn.nstateful
    ind = sn.statefulToNode(isf);
    if size(state{isf},1) > 1
        line_warning(mfilename,'Some states at node %d will be ignored. Please assign the node with a specific state.\n', ind);
    end
    [~,nir(isf,:)] = State.toMarginal(sn, ind, state{isf}(1,:));
end
nir = nir';
rows = findrows(TSS(:,2:end), nir(:)');
if ~isempty(rows)
    ProbSysAggr = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.\n');
    ProbSysAggr = 0;
end
end