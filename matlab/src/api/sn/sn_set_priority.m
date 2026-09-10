%{ @file sn_set_priority.m
 %  @brief Sets the priority for a class
 %
 %  @author LINE Development Team
%}

%{
 % @brief Sets the priority for a class
 %
 % @details
 % Directly modifies the class priority in NetworkStruct.
 %
 % @par Syntax:
 % @code
 % sn = sn_set_priority(sn, classIdx, priority)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>classIdx<td>Class index (1-based)
 % <tr><td>priority<td>Priority value
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Modified network structure
 % </table>
%}
function sn = sn_set_priority(sn, classIdx, priority)

sn.classprio(classIdx) = priority;

end
