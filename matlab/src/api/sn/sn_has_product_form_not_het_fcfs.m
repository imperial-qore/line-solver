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
function bool = sn_has_product_form_not_het_fcfs(sn)

bool = all(sn.sched==SchedStrategy.INF | sn.sched==SchedStrategy.PS | sn.sched==SchedStrategy.FCFS | sn.sched==SchedStrategy.LCFSPR | sn.sched==SchedStrategy.EXT);
bool = bool & ~sn_has_priorities(sn);
bool = bool & ~sn_has_fork_join(sn);
bool = bool & ~sn_has_sd_routing(sn);
iset = find(sn.sched == SchedStrategy.FCFS);
for i=iset(:)'
    icset = isfinite(sn.scv(i,:)) & sn.scv(i,:)>0;
    bool = bool & all(sn.scv(i,icset) > 1-GlobalConstants.FineTol) & all(sn.scv(i,icset) < 1+GlobalConstants.FineTol);
end
end