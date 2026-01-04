%{ @file sn_has_multi_class_heter_fcfs.m
 %  @brief Checks for multi-class FCFS with heterogeneous service
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks for multi-class FCFS with heterogeneous service
 %
 % @details
 % Returns true if any FCFS station serves multiple classes with different service time distributions.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_multi_class_heter_fcfs(sn)
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
 % <tr><td>bool<td>True if any multi-class FCFS station has heterogeneous service
 % </table>
%}
function bool = sn_has_multi_class_heter_fcfs(sn)

iset = find(sn.sched == SchedStrategy.FCFS);
if isempty(iset)
    bool = false;
else
    bool = false;
    for i=iset(:)'
        bool = bool | range([sn.rates(i,:)])>0;
    end
end
end
