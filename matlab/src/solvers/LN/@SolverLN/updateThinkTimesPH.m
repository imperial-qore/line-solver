function updateThinkTimesPH(self, it)
% UPDATETHINKTIMESPH(IT) Surrogate delay of every caller under method 'srvn.ph'
%
% Same closure as updateThinkTimes -- a thread of the task is idle for whatever
% of its cycle the task's own station does not hold -- but the rate it is
% normalised by is the INVOCATION rate of the task and not the throughput of its
% station. Under this method a caller class reaches the server once per
% invocation of the caller, carrying its whole call burst in its service law, so
% the station rate counts caller cycles rather than calls and the two differ by
% the mean number of calls.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

lqn = self.lqn;

for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if self.ignore(tidx)
        continue
    end
    % only a reference task's think time separates one request from the next;
    % on a served task it is not a per-request delay -- see lqn_ref_thinktime
    ztask = lqn_ref_thinktime(lqn, tidx);
    if isnan(self.idxhash(tidx))
        % A task no other task calls but whose entries carry an arrival still has
        % a cycle: its threads are driven by the stream. buildLayersPH drops the
        % open class for it precisely so this closure can set the rate, the way
        % updateThinkTimes does for the default method. Without it the chain runs
        % against an Immediate delay and saturates -- lqn_open_arrival reported
        % the processor at 1.000 where lqns, lqsim and LDES all give 0.32.
        arvrate = phArrivalRate(lqn, tidx);
        if arvrate > GlobalConstants.FineTol
            njobs = lqn.maxmult(tidx);
            if ~isfinite(njobs) || njobs <= 0
                njobs = max(self.njobs(tidx,:));
            end
            z = max(GlobalConstants.Zero, njobs/arvrate - phHostResid(self, lqn, tidx) - ztask);
            omega = self.relax_omega;
            if omega < 1.0 && it > 1 && ~isnan(self.thinkt_prev(tidx))
                z = omega * z + (1 - omega) * self.thinkt_prev(tidx);
            end
            self.tput(tidx) = arvrate;
            self.thinkt(tidx) = z;
            self.thinkt_prev(tidx) = z;
            self.thinktproc{tidx} = Exp.fitMean(z + ztask);
            continue
        end
        % a reference task, or one no other task calls: it has no station of
        % its own, so its only delay is the think time the user declared
        self.thinkt(tidx) = GlobalConstants.FineTol;
        self.thinktproc{tidx} = Immediate();
        continue
    end
    [layer_t, station_t] = self.layerOf(tidx);
    U = sum(self.results{end,layer_t}.UN(station_t,:), 2);
    self.util(tidx) = U;
    % The closure below measures a thread's cycle in WORK: it is idle for
    % whatever of the cycle its station does not hold it working. A SetupTask's
    % station service also carries a cold start, which is time the thread is
    % unavailable but is not work, so it is taken back out of U before the
    % closure reads it. Leaving it in shortened the surrogate delay and ran the
    % callee layer 5.8% above the rate its callers actually drive it at, on
    % lqn_setup. Zero for every task without a setup.
    if isfield(self.ph,'setupshare') && ~isempty(self.ph.setupshare) ...
            && tidx <= numel(self.ph.setupshare) && self.ph.setupshare(tidx) > 0
        U = U * (1 - self.ph.setupshare(tidx));
    end
    % the rate the CALLERS ask of the task, not the rate its processor layer
    % reported: the latter is itself a function of this think time
    X = self.ph.xdemand(tidx);
    if ~(X > GlobalConstants.FineTol)
        X = self.tput(tidx);
    end
    % The thread pool of ONE replica, the convention ph.xdemand is kept in.
    % Reading it off the layers takes whichever convention each layer happened
    % to use for this caller: on sockshop the task layer of T1's callee counts
    % 48 threads, both replicas, against a per-replica rate, and the closure
    % then reports T1 at half the rate its own caller drives it at.
    njobs = lqn.maxmult(tidx);
    if ~isfinite(njobs) || njobs <= 0
        njobs = max(self.njobs(tidx,:));
    end
    if X > GlobalConstants.FineTol
        if lqn.sched(tidx) == SchedStrategy.INF
            % an infinite server reports a mean number of busy threads
            z = (njobs - U) / X - ztask;
        else
            z = njobs * abs(1 - U) / X - ztask;
        end
    else
        z = self.thinkt(tidx);
    end
    z = max(GlobalConstants.Zero, z);
    if it > 1 && ~isnan(self.thinkt_prev(tidx)) && (isinf(z) || isnan(z))
        z = self.thinkt_prev(tidx);
    end
    omega = self.relax_omega;
    if omega < 1.0 && it > 1 && ~isnan(self.thinkt_prev(tidx))
        z = omega * z + (1 - omega) * self.thinkt_prev(tidx);
    end
    self.thinkt(tidx) = z;
    self.thinkt_prev(tidx) = z;
    self.thinktproc{tidx} = Exp.fitMean(z + ztask);
end
end

% ------------------------------------------------------------------------
function rate = phArrivalRate(lqn, tidx)
% Total exogenous rate into the entries of task TIDX, zero unless the arrival is
% the only way in -- the predicate buildLayersPH drops the open class on
rate = 0;
if lqn.isref(tidx)
    return
end
for eidx = lqn.entriesof{tidx}
    if any(full(lqn.issynccaller(:, eidx))) || any(full(lqn.isasynccaller(:, eidx)))
        return
    end
end
if ~(isfield(lqn,'arrival') && ~isempty(lqn.arrival) && iscell(lqn.arrival))
    return
end
for eidx = lqn.entriesof{tidx}
    if eidx <= numel(lqn.arrival) && ~isempty(lqn.arrival{eidx})
        m = lqn.arrival{eidx}.getMean();
        if isfinite(m) && m > GlobalConstants.FineTol
            rate = rate + 1/m;
        end
    end
end
end

% ------------------------------------------------------------------------
function r = phHostResid(self, lqn, tidx)
% Response time the caller class of task TIDX sees at its processor layer
r = 0;
hidx = lqn.parent(tidx);
if isnan(hidx) || hidx < 1 || hidx > numel(self.idxhash) || isnan(self.idxhash(hidx))
    return
end
lay = self.idxhash(hidx);
tasks = self.ensemble{lay}.attribute.tasks;
row = find(tasks(:,2) == tidx & ~isnan(tasks(:,1)), 1);
if isempty(row)
    return
end
stn = self.ensemble{lay}.attribute.serverIdxOf(hidx);
if isnan(stn)
    return
end
res = self.results{end, lay};
if isempty(res) || isnan(res.RN(stn, tasks(row,1)))
    return
end
r = res.RN(stn, tasks(row,1));
end
