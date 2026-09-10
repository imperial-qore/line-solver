%{ @file sn_is_population_model.m
 %  @brief Checks if the network is defined by populations
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network is defined by populations
 %
 % @details
 % Returns true if the model is specified using population values rather than arrival rates.
 %
 % @par Syntax:
 % @code
 % bool = sn_is_population_model(sn)
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
 % <tr><td>bool<td>True if model is population-based
 % </table>
%}
function bool = sn_is_population_model(sn)

bool = all(sn.sched==SchedStrategy.INF | sn.sched==SchedStrategy.PS ...
    | sn.sched==SchedStrategy.PSPRIO  | sn.sched==SchedStrategy.DPS ...
    | sn.sched==SchedStrategy.GPS  | sn.sched==SchedStrategy.GPSPRIO ...
    | sn.sched==SchedStrategy.DPSPRIO | sn.sched==SchedStrategy.EXT);
bool = bool & ~sn_has_priorities(sn);
bool = bool & ~sn_has_fork_join(sn);
end