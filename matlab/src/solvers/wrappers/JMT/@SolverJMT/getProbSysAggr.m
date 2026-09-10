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

function probSysStateAggr = getProbSysAggr_impl(self)
% GETPROBSYSAGGR_IMPL Implementation of GETPROBSYSAGGR; see the wrapper above.
if GlobalConstants.DummyMode
    probSysStateAggr = NaN;
    return
end
% lang='cpp' answers this from `-s jmt -a prob`; see getProbAggr for why the
% delegation is the same estimator and still the C++ engine's number.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    CPPLINE.assertSingleState(self.name, 'getProbSysAggr', self.model);
    payload = CPPLINE.analysisViaCpp(self.name, self.model, self.options, 'prob');
    if ~isfield(payload, 'ProbSysAggr')
        CPPLINE.cppUnsupported(self.name, 'getProbSysAggr', ...
            'line-cli''s -a prob answer for this solver carries no ProbSysAggr field');
    end
    probSysStateAggr = CPPLINE.jsonNumericScalar(payload.ProbSysAggr);
    % The reference warns when the state never occurred, because on a simulation
    % a zero is far more often the run length than the model; line-cli reports
    % that fact beside the number, so the warning can be the same one.
    if isfield(payload, 'SysStateSeen') && ~payload.SysStateSeen
        line_warning(mfilename,'The state was not seen during the simulation.\n');
    end
    return
end

sn = self.getStruct;
TranSysStateAggr = self.sampleSysAggr;
TSS = cell2mat([TranSysStateAggr.t,TranSysStateAggr.state(:)']);
TSS(:,1)=[diff(TSS(:,1));0];
state = sn.state;
nir = zeros(sn.nstateful,sn.nclasses);
for isf=1:sn.nstateful
    ind = sn.statefulToNode(isf);
    [~,nir(isf,:)] = State.toMarginal(sn, ind, state{isf});
end
nir = nir';
rows = findrows(TSS(:,2:end), nir(:)');
if ~isempty(rows)
    probSysStateAggr = sum(TSS(rows,1))/sum(TSS(:,1));
else
    line_warning(mfilename,'The state was not seen during the simulation.\n');
    probSysStateAggr = 0;
end
end