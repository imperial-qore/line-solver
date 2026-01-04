function CDc = getCdfPT(self)
% CDC = GETCDFPT()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0 = tic;
%if ~exist('R','var')
R = self.getAvgRespTHandles;
%end
sn = self.getStruct;
%if ~exist('withRefStat','var')
withRefStat = false(1,sn.nclasses);
%elseif numel(withRefStat) == 1
%    val = withRefStat;
%    withRefStat = false(1,sn.nclasses);
%    withRefStat(:) = val;
%end
% ptSpec = struct(); % passage time specification
% ptSpec.starts = false(sn.nnodes,sn.nclasses,sn.nnodes,sn.nclasses);
% ptSpec.completes = false(sn.nnodes,sn.nclasses,sn.nnodes,sn.nclasses);
% for r=1:sn.nclasses
%     if withRefStat(r)
%         % starts when arriving to ref
%         ptSpec.starts(:,:,sn.refstat(r),r) = true;
%     else % ref station excluded
%         % starts when leaving ref
%         ptSpec.starts(sn.refstat(r),r,:,:) = true;
%     end
%     % completes when arriving to ref
%     if R{sn.refstat(r),r}.class.completes
%         % class switch to r is right after departure from station i
%         ptSpec.completes(:,:,sn.refstat(r),r) = true;
%     end
% end
% options = self.getOptions;
% options.psgtime = ptSpec;
completes = false(sn.nnodes,sn.nclasses);
for i=1:sn.nstations
    for r=1:sn.nclasses
        if R{i,r}.class.completes
            completes(i,r) = true;
        end
    end
end
CDc = solver_fluid_RT(sn, self.result.solverSpecific.odeStateVec, options, completes);
runtime = toc(T0);
self.setDistribResults(CDc, runtime);
end
