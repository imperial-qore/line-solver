%{ @file sn_has_bursty_arrival.m
 %  @brief Checks whether the network has a bursty (non-renewal) arrival process
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks whether any external arrival process is bursty, i.e. non-renewal
 %
 % @details
 % Returns true if any Source station has an arrival process with autocorrelated
 % inter-arrival times (a non-renewal Markovian arrival process such as an
 % MMPP/MAP), as opposed to a renewal process (Poisson, or any i.i.d. renewal
 % process such as Erlang/HyperExp/Coxian/APH). Detection is exact: a MAP with
 % matrices (D0,D1) is renewal iff D1 equals its rank-one renewal form t0*pie,
 % where t0 = -D0*e and pie is the embedded stationary vector; any departure from
 % that form signals correlation between successive inter-arrival times.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_bursty_arrival(sn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>bool<td>True if some external arrival process is non-renewal (bursty)
 % </table>
%}
function bool = sn_has_bursty_arrival(sn)

bool = false;
for ist = 1:sn.nstations
    nd = sn.stationToNode(ist);
    if sn.nodetype(nd) ~= NodeType.Source
        continue;
    end
    for r = 1:sn.nclasses
        % sn.proc is a per-station cell of per-class representations, so the
        % class index is the SECOND brace: sn.proc{ist,r} reads the station cell
        % and only ever reached class 1 (and threw for r>1 on a multiclass
        % model, which the nclasses==1 gates of its callers hid).
        if isempty(sn.proc) || numel(sn.proc) < ist || numel(sn.proc{ist}) < r
            continue;
        end
        map = sn.proc{ist}{r};
        if isempty(map) || numel(map) < 2 || isempty(map{1})
            continue;
        end
        D1 = map{2};
        n = size(D1,1);
        if n <= 1
            continue;   % single-phase arrival is Poisson, hence renewal
        end
        D1ren = D1 * ones(n,1) * map_pie(map);   % rank-one renewal form t0*pie
        if norm(D1 - D1ren, 'fro') > 1e-8 * max(1, norm(D1,'fro'))
            bool = true;
            return;
        end
    end
end
end
