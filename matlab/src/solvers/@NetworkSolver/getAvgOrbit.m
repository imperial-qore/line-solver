function ON = getAvgOrbit(self)
% ON = GETAVGORBIT()
%
% Mean number of jobs waiting in the ORBIT of each retrial station, as an
% (nstations x nclasses) matrix. Stations that are not retrial queues report 0.
%
% A retrial station has no waiting room: a job that finds every server busy
% joins the orbit instead of queueing, so its station population splits into
% the jobs currently in service and the jobs orbiting. getAvgQLen reports the
% whole station population, which is why the orbit had to be recovered by hand
% as QLen - Util. This method reports it directly.
%
% The in-service population is obtained from the station throughput by Little's
% law applied to the servers alone, E[in service] = X * E[S], which holds for
% any service distribution and any number of servers, so the orbit length is
% exact whenever QLen and Tput are.
%
% See also Queue.setOrbit, NetworkSolver.getAvgTable
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[QN,~,~,TN] = self.getAvg();
sn = self.model.getStruct();
ON = zeros(size(QN));

if ~isfield(sn,'retrialProc') || isempty(sn.retrialProc)
    return
end

for ist = 1:size(ON,1)
    for r = 1:size(ON,2)
        if size(sn.retrialProc,1) < ist || size(sn.retrialProc,2) < r ...
                || isempty(sn.retrialProc{ist,r})
            continue % not a retrial station for this class: no orbit
        end
        rate_ir = sn.rates(ist,r);
        if ~isfinite(rate_ir) || rate_ir <= 0
            continue
        end
        inService = TN(ist,r) / rate_ir; % Little's law on the servers
        ON(ist,r) = max(0, QN(ist,r) - inService);
    end
end
end
