function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)

T0 = tic;
if nargin<2 %~exist('R','var')
    R = self.getAvgRespTHandles;
end
sn = self.getStruct;
self.getAvg; % get steady-state solution
options = self.getOptions;
RD = solver_mam_passage_time(sn, sn.proc, options);
runtime = toc(T0);
self.setDistribResults(RD, runtime);
end