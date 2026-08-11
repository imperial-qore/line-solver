%{ @file sn_has_product_form.m
 %  @brief Checks if the network has product-form solution
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has product-form solution
 %
 % @details
 % Returns true if the network satisfies conditions for product-form equilibrium distribution.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_product_form(sn)
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
 % <tr><td>bool<td>True if the network has product-form solution
 % </table>
%}
function bool = sn_has_product_form(sn)

bool = all(sn.sched==SchedStrategy.INF | sn.sched==SchedStrategy.PS | sn.sched==SchedStrategy.FCFS | sn.sched==SchedStrategy.LCFSPR | sn.sched==SchedStrategy.LCFS | sn.sched==SchedStrategy.EXT);
bool = bool & ~sn_has_multi_class_heter_fcfs(sn);
bool = bool & ~sn_has_priorities(sn);
bool = bool & ~sn_has_fork_join(sn);
bool = bool & ~sn_has_sd_routing(sn);
end