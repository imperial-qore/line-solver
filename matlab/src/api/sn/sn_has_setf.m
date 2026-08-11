%{ @file sn_has_setf.m
 %  @brief Checks if the network uses SETF scheduling
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network uses SETF scheduling
 %
 % @details
 % Returns true if any station in the network uses Shortest Elapsed Time First
 % (SETF) scheduling. SETF is the non-preemptive version of FB/LAS scheduling.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_setf(sn)
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
 % <tr><td>bool<td>True if any station uses SETF scheduling
 % </table>
%}
function bool = sn_has_setf(sn)

bool = any(sn.sched==SchedStrategy.SETF);
end
