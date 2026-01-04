%{ @file sn_set_fork_fanout.m
 %  @brief Sets fork fanout (tasksPerLink) for a Fork node
 %
 %  @author LINE Development Team
%}

%{
 % @brief Sets fork fanout for a Fork node
 %
 % @details
 % Updates the fanOut field in nodeparam for a Fork node.
 %
 % @par Syntax:
 % @code
 % sn = sn_set_fork_fanout(sn, forkNodeIdx, fanOut)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>forkNodeIdx<td>Node index of the Fork node (1-based)
 % <tr><td>fanOut<td>Number of tasks per output link (>= 1)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Modified network structure
 % </table>
%}
function sn = sn_set_fork_fanout(sn, forkNodeIdx, fanOut)

% Verify it's a Fork node
if sn.nodetype(forkNodeIdx) ~= NodeType.Fork
    error('sn_set_fork_fanout: Node %d is not a Fork node', forkNodeIdx);
end

% Update nodeparam
if ~isfield(sn.nodeparam{forkNodeIdx}, 'fanOut')
    sn.nodeparam{forkNodeIdx}.fanOut = fanOut;
else
    sn.nodeparam{forkNodeIdx}.fanOut = fanOut;
end

end
