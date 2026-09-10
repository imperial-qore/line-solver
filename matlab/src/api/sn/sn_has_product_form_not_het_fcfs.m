%{ @file sn_has_product_form_not_het_fcfs.m
 %  @brief Checks for product-form excluding heterogeneous FCFS
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks for product-form excluding heterogeneous FCFS
 %
 % @details
 % Returns true if the network has product-form solution excluding cases with heterogeneous FCFS.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_product_form_not_het_fcfs(sn)
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
 % <tr><td>bool<td>True if product-form holds without heterogeneous FCFS
 % </table>
%}
function bool = sn_has_product_form_not_het_fcfs(sn, checkMeans)
% CHECKMEANS (default true) also demands class-independent FCFS service means.
% Pass false only for an algorithm that models class-dependent FCFS itself
% (ab, schmidt, schmidt-ext), for which the exclusion is the whole point.
if nargin < 2
    checkMeans = true;
end

bool = all(sn.sched==SchedStrategy.INF | sn.sched==SchedStrategy.PS | sn.sched==SchedStrategy.FCFS | sn.sched==SchedStrategy.LCFSPR | sn.sched==SchedStrategy.EXT);
bool = bool & ~sn_has_priorities(sn);
bool = bool & ~sn_has_fork_join(sn);
bool = bool & ~sn_has_sd_routing(sn);
iset = find(sn.sched == SchedStrategy.FCFS);
for i=iset(:)'
    icset = isfinite(sn.scv(i,:)) & sn.scv(i,:)>0;
    bool = bool & all(sn.scv(i,icset) > 1-GlobalConstants.FineTol) & all(sn.scv(i,icset) < 1+GlobalConstants.FineTol);
    % BCMP type 1 asks the FCFS service to be exponential AND class-independent,
    % so the means must agree too: with unequal means the product-form solve
    % returns a wait proportional to each class's own demand where FCFS makes
    % every class wait behind the same queue. The comparison is between CHAIN
    % service times (visit-weighted over the classes that actually visit the
    % station): a class that never visits cannot break product form, and
    % within-chain heterogeneity is invisible to both the product-form and the
    % qd branch, which deaggregate a chain result proportionally to each
    % class's own demand, so only between-chain heterogeneity warrants the
    % divert. LN layers carry seeded rates for classes with zero visits, which
    % a raw per-class comparison mistakes for heterogeneity.
    if checkMeans
        isf = sn.stationToStateful(i);
        stchain = [];
        for c = 1:sn.nchains
            inchain = find(sn.chains(c,:));
            w = full(sn.visits{c}(isf,inchain));
            imset = w > GlobalConstants.Zero & isfinite(sn.rates(i,inchain)) & sn.rates(i,inchain) > 0;
            if any(imset)
                stchain(end+1) = sum(w(imset) ./ sn.rates(i,inchain(imset))) / sum(w(imset)); %#ok<AGROW>
            end
        end
        if ~isempty(stchain)
            bool = bool & (max(stchain) - min(stchain) <= GlobalConstants.CoarseTol * max(stchain));
        end
    end
end
end