%{
 % @brief Majumdar-Woodside robust box bounds on throughput for a layered
 %        queueing network (LQN), computed on the processor-contention model.
 %
 % @details
 % Builds a closed multiclass queueing network from a LayeredNetworkStruct in
 % which the stations are the processors (hosts) and the classes are the
 % reference-task call chains, then applies the Majumdar-Woodside robust box
 % bounds (see pfqn_mwrbb). The per-chain demand at each processor is the total
 % host demand executed on that processor during one cycle of the reference
 % task, obtained by traversing the activity/call graph and scaling by the
 % mean number of synchronous calls. The reference-task multiplicity is the
 % class population and the reference-task think time is the class think time.
 %
 % This generalizes the classical LQN "Type 1 throughput bound"
 % X_ref <= mult/(Z + D_total) (the no-contention upper bound, as computed by
 % lqns -b) by additionally providing the processor-utilization upper bound and
 % the Majumdar-Woodside lower bound (throughput guarantee, Theorem 2).
 %
 % Scope: processors are the queueing resources; reference-task chains are the
 % classes. Software (finite-thread task) bottlenecks and non-deterministic
 % activity precedence (OR-branch probabilities, loop repetitions) are not
 % modeled here (sequential activity execution is assumed).
 %
 % @par Syntax:
 % @code
 % out = lqn_boxbounds(lqn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lqn<td>LayeredNetworkStruct (from LayeredNetwork.getStruct())
 % </table>
 %
 % @par Returns (struct out):
 % <table>
 % <tr><th>Field<th>Description
 % <tr><td>refidx<td>(1 x R) absolute indices of the reference tasks
 % <tr><td>Xlo,Xup<td>(1 x R) lower/upper throughput bound per reference chain
 % <tr><td>TN_lo,TN_up<td>(nidx x 1) throughput bound propagated to all elements
 % <tr><td>UN_lo,UN_up<td>(nidx x 1) processor utilization bound (util law)
 % <tr><td>D<td>(nhosts x R) per-chain demand at each processor
 % </table>
%}
function out = lqn_boxbounds(lqn)
nidx = lqn.nidx;
nH = lqn.nhosts;

% reference tasks -> classes
refidx = [];
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if lqn.isref(tidx)
        refidx(end+1) = tidx; %#ok<AGROW>
    end
end
R = numel(refidx);

D = zeros(nH, R);
Vis = zeros(nidx, R);
Nref = zeros(1, R);
Zref = zeros(1, R);
for r = 1:R
    tidx = refidx(r);
    Nref(r) = lqn.mult(tidx);
    if isinf(Nref(r)), Nref(r) = 1; end     % open/infinite ref treated as 1
    z = lqn.think_mean(tidx);
    if isnan(z), z = 0; end
    Zref(r) = z;
    for eidx = lqn.entriesof{tidx}
        [D(:,r), Vis(:,r)] = visitEntry(lqn, eidx, 1.0, nH, D(:,r), Vis(:,r));
    end
end

% processor scheduling -> discipline code; equal class priority
schedH = zeros(nH, 1);
for h = 1:nH
    schedH(h) = discCode(lqn.sched(h));
end
prio = zeros(1, R);

V = double(D > 0);
S = D;
[Xlo, Xup] = pfqn_mwrbb(V, S, Nref, Zref, schedH, prio);

% propagate to all LQN elements
TN_lo = nan(nidx, 1); TN_up = nan(nidx, 1);
UN_lo = nan(nidx, 1); UN_up = nan(nidx, 1);
for idx = 1:nidx
    if any(Vis(idx, :) > 0)
        TN_lo(idx) = sum(Xlo .* Vis(idx, :));
        TN_up(idx) = sum(Xup .* Vis(idx, :));
    end
end
for h = 1:nH
    UN_lo(h) = sum(Xlo .* D(h, :));
    UN_up(h) = sum(Xup .* D(h, :));
end

out = struct('refidx', refidx, 'Xlo', Xlo, 'Xup', Xup, ...
    'TN_lo', TN_lo, 'TN_up', TN_up, 'UN_lo', UN_lo, 'UN_up', UN_up, 'D', D);
end

% ------------------------------------------------------------------------
function [d, vis] = visitEntry(lqn, eidx, mult, nH, d, vis)
vis(eidx) = vis(eidx) + mult;
tidx = lqn.parent(eidx);           % task hosting this entry
if tidx >= 1 && tidx <= numel(vis)
    vis(tidx) = vis(tidx) + mult;
end
acts = lqn.actsof{eidx};
for aidx = acts
    if lqn.parent(aidx) == lqn.parent(eidx)   % activity of this entry
        [d, vis] = visitActivity(lqn, aidx, mult, nH, d, vis);
    end
end
end

function [d, vis] = visitActivity(lqn, aidx, mult, nH, d, vis)
vis(aidx) = vis(aidx) + mult;
tidx = lqn.parent(aidx);
hidx = lqn.parent(tidx);              % host absolute index (hshift=0)
hd = lqn.hostdem_mean(aidx);
if isnan(hd), hd = 0; end
if hidx >= 1 && hidx <= nH
    d(hidx) = d(hidx) + mult * hd;
end
for cidx = lqn.callsof{aidx}
    if lqn.calltype(cidx) == CallType.SYNC
        cmean = lqn.callproc_mean(cidx);
        callee = lqn.callpair(cidx, 2);
        [d, vis] = visitEntry(lqn, callee, mult * cmean, nH, d, vis);
    end
end
end

% Map a SchedStrategy id to a Majumdar-Woodside discipline code:
% 0=FIFO, 1=PS, 2=non-preemptive priority, 3=preemptive priority,
% 4=ABA full-contention (discipline-independent).
function code = discCode(s)
switch s
    case {SchedStrategy.FCFS}
        code = 0;
    case {SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, ...
          SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO}
        code = 1;
    case {SchedStrategy.HOL}
        code = 2;
    case {SchedStrategy.FCFSPRPRIO, SchedStrategy.LCFSPRPRIO}
        code = 3;
    otherwise
        code = 4;
end
end
