function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)

sn = self.getStruct;
if GlobalConstants.DummyMode
    RD = cell(sn.nstations,sn.nclasses);
    logData = NaN;
    return
end

T0 = tic;
if nargin<2 %~exist('R','var')
    R = self.getAvgRespTHandles;
    % to do: check if some R are disabled
end
self.getAvg; % get steady-state solution
options = self.getOptions;
options.init_sol = self.result.solverSpecific.odeStateVec;
% we need to pass the modified sn as the number of phases may have changed
% during the fluid iterations, affecting the size of odeStateVec
RD = solver_fluid_passage_time(self.result.solverSpecific.sn, options);
runtime = toc(T0);
self.setDistribResults(RD, runtime);
end