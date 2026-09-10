function RD = getCdfPassT(self, R)
% RD = GETCDFPASST(R)
%
% Passage time distribution. For SolverMAM the passage time IS the response
% time: both come from the same solver_mam_passage_time call, so this
% delegates to getCdfRespT, as the JAR and the C++ CLI arms do.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    RD = self.getCdfRespT();
else
    RD = self.getCdfRespT(R);
end
end
