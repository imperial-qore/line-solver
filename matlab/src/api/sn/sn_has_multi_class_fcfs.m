%{ @file sn_has_multi_class_fcfs.m
 %  @brief Checks if the network has multi-class FCFS stations
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has multi-class FCFS stations
 %
 % @details
 % Returns true if any FCFS station serves multiple job classes.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_multi_class_fcfs(sn)
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
 % <tr><td>bool<td>True if any FCFS station is multi-class
 % </table>
%}
function bool = sn_has_multi_class_fcfs(sn)

iset = find(sn.sched == SchedStrategy.FCFS);
bool = false;
for i=iset(:)'
    bool = bool | sum(sn.rates(i,:)>0)>1;
end
end