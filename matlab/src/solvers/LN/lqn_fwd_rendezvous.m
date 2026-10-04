function lqn = lqn_fwd_rendezvous(lqn)
% LQN = LQN_FWD_RENDEZVOUS(LQN)
% Forwarding transformation of Franks (1999), Sec. 3.3.1 and Fig. 3.8: each
% forwarding chain reachable from a synchronous call is reconnected to the
% client that issued the original rendezvous, as a pseudo rendezvous (SYNC)
% call whose mean is the original call mean times the product of the
% forwarding probabilities along the path. One level of servers disappears
% from the layering, and the forwarded workload is then carried by ordinary
% SYNC call classes, so layer construction, think times, populations and the
% interlock analysis all see plain rendezvous arcs.
%
% As the thesis notes, the transformed model is not one in which the client
% makes two remote procedure calls directly: the pseudo arcs are excluded
% from the slice times and from the overtaking and interlock probabilities.
% FWD calls are kept in the struct but no longer contribute blocking
% anywhere in SolverLN.
%
% Asynchronous calls into a forwarding chain are left untouched, since a
% send-no-reply terminates the chain of blocking.

if ~any(lqn.calltype == CallType.FWD)
    return
end

ncalls0 = lqn.ncalls;
for cidx = 1:ncalls0
    if lqn.calltype(cidx) ~= CallType.SYNC
        continue
    end
    aidx = lqn.callpair(cidx,1);
    tidx = lqn.parent(aidx);
    base_mean = lqn.callproc_mean(cidx);
    if base_mean <= 0
        continue
    end
    % BFS through the forwarding chain of the sync target
    frontier = lqn.callpair(cidx,2);
    probs = 1;
    visitedE = [];
    while ~isempty(frontier)
        eidx = frontier(1); frontier(1) = [];
        p_path = probs(1); probs(1) = [];
        if ismember(eidx, visitedE), continue; end
        visitedE(end+1) = eidx; %#ok<AGROW>
        for fcidx = 1:ncalls0
            if lqn.calltype(fcidx) ~= CallType.FWD || lqn.callpair(fcidx,1) ~= eidx
                continue
            end
            fprob = lqn.callproc_mean(fcidx);
            tgt = lqn.callpair(fcidx,2);
            pseudo_mean = base_mean * p_path * fprob;
            if pseudo_mean > 0 && lqn.parent(tgt) ~= tidx
                % see _kb/06-solver-catalog.md (LN section) for rationale
                mrow = 0;
                for scan = 1:lqn.ncalls
                    if lqn.calltype(scan) == CallType.SYNC && ...
                            lqn.callpair(scan,1) == aidx && lqn.callpair(scan,2) == tgt
                        mrow = scan;
                        break
                    end
                end
                if mrow > 0
                    newmean = lqn.callproc_mean(mrow) + pseudo_mean;
                    d = Geometric(1/newmean);
                    lqn.callproc{mrow,1} = d;
                    lqn.callproc_mean(mrow) = newmean;
                    lqn.callproc_scv(mrow) = d.getSCV();
                else
                    ncall = lqn.ncalls + 1;
                    lqn.ncalls = ncall;
                    target_tidx = lqn.parent(tgt);
                    d = Geometric(1/pseudo_mean);
                    lqn.calltype(ncall,1) = CallType.SYNC;
                    lqn.callpair(ncall,1:2) = [aidx, tgt];
                    lqn.callnames{ncall,1} = [lqn.names{aidx},'=>',lqn.names{tgt}];
                    lqn.callhashnames{ncall,1} = [lqn.hashnames{aidx},'=>',lqn.hashnames{tgt}];
                    lqn.callproc{ncall,1} = d;
                    % Only the mean (and the process object) are consumed by
                    % SolverLN; mirror the base call for the remaining fields
                    lqn.callproc_type(ncall) = lqn.callproc_type(cidx);
                    lqn.callproc_params{ncall} = lqn.callproc_params{cidx};
                    lqn.callproc_mean(ncall) = pseudo_mean;
                    lqn.callproc_scv(ncall) = d.getSCV();
                    lqn.callproc_proc{ncall} = lqn.callproc_proc{cidx};
                    lqn.callsof{aidx}(end+1) = ncall;
                    lqn.iscaller(tidx, target_tidx) = true;
                    lqn.iscaller(aidx, target_tidx) = true;
                    lqn.iscaller(tidx, tgt) = true;
                    lqn.iscaller(aidx, tgt) = true;
                    lqn.issynccaller(tidx, target_tidx) = true;
                    lqn.issynccaller(aidx, target_tidx) = true;
                    lqn.issynccaller(tidx, tgt) = true;
                    lqn.issynccaller(aidx, tgt) = true;
                    lqn.graph(aidx, tgt) = 1;
                    lqn.taskgraph(tidx, target_tidx) = 1;
                end
            end
            % Follow the chain
            if ~ismember(tgt, visitedE) && ~ismember(tgt, frontier)
                frontier(end+1) = tgt; %#ok<AGROW>
                probs(end+1) = p_path * fprob; %#ok<AGROW>
            end
        end
    end
end
end
