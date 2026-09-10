function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)


% lang='cpp' takes the law from line-cli (-s mam -a cdf), which runs the same
% MMAPPH1-family passage-time analysis solver_mam_passage_time runs below.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    RD = CPPLINE.mamCdfRespT(self.name, self.model, self.options);
    return
end

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