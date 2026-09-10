function varargout = getProbSys(self,varargin)
% PROBSYS = GETPROBSYS()
% The result recorder captures the scalar this getter returned together
% with the solver that produced it -- see LineResultRecorder. Six of the
% statepr_* goldens hold exactly this number and nothing else, so recording
% it is what makes those goldens attributable instead of "the first bare
% number the example printed". The wrapper exists so that recording happens
% on EVERY exit path of the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getProbSys_impl(self,varargin{:});
LineResultRecorder.captureScalar(scope, self, 'probSys', varargout{1});
end

function probSys = getProbSys_impl(self)
% GETPROBSYS_IMPL Implementation of GETPROBSYS; see the wrapper above.

% lang='cpp' answers this from -a prob. See CPPLINE.assertSingleState.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    CPPLINE.assertSingleState(self.name, 'getProbSys', self.model);
    probSys = CPPLINE.probAggr(self.name, self.model, self.options).ProbSys;
    return
end

if GlobalConstants.DummyMode
    probSys = NaN;
    return
end
TranSysState = self.sampleSys;
TSS = cell2mat([TranSysState.t,TranSysState.state(:)']);
TSS(:,1)=[TSS(1,1);diff(TSS(:,1))];
sn = self.getStruct;
% fill-in for FCFS states
for isf=1:size(TranSysState.state,2)
    if size(sn.state{isf},1)>1
        error('There are multiple station states, choose an initial state as a parameter to getProb.');
    end
    sn.state{isf} = [zeros(1,size(TranSysState.state{isf},2)-size(sn.state{isf},2)),sn.state{isf}];
end
state = cell2mat(sn.state');
rows = findrows(TSS(:,2:end), state);
if ~isempty(rows)
    probSys = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.\n');
    probSys = 0;
end
end