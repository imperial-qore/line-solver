%{ @file sn_has_multi_class_heter_exp_fcfs.m
 %  @brief Checks for multi-class FCFS with heterogeneous exponential service
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks for multi-class FCFS with heterogeneous exponential service
 %
 % @details
 % Returns true if any FCFS station serves multiple classes with different exponential service times.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_multi_class_heter_exp_fcfs(sn)
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
 % <tr><td>bool<td>True if any multi-class FCFS station has heterogeneous exponential service
 % </table>
%}
function bool = sn_has_multi_class_heter_exp_fcfs(sn)
% Checks if the network has one or more stations with multiclass heterogeneous FCFS
% and exponential service times
%
% Parameters:
%   sn - NetworkStruct object for the queueing network model
%
% Returns:
%   bool - true if network has multiclass heterogeneous FCFS with exponential service

bool = false;
iset = find(sn.sched == SchedStrategy.FCFS);

if isempty(iset)
    return;
end

for i = iset(:)'
    % Check if rates vary across classes (heterogeneous)
    if range(sn.rates(i,:)) > 0
        % Check if all SCVs are ~1 (exponential)
        scvs = sn.scv(i,:);
        if max(scvs) < 1 + GlobalConstants.FineTol && ...
           min(scvs) > 1 - GlobalConstants.FineTol
            bool = true;
            return;
        end
    end
end
end