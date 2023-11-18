function refreshJobs(self)
% REFRESHJOBS()
% Updates sn structure after a change in a closed class population

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

refreshStruct(self);
njobs = getNumberOfJobs(self);
self.sn.nclosedjobs = sum(njobs(isfinite(njobs)));
self.sn.njobs = njobs(:)';
end
