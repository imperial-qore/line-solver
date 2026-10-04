%{ @file sn_has_quorum_join.m
 %  @brief Checks if the network has a quorum (k-of-n) join
 %
 %  @author LINE Development Team
%}

%{
 % @brief Checks if the network has a quorum (k-of-n) join
 %
 % @details
 % Returns true if some Join node declares a non-standard strategy with a
 % positive required count in some class, i.e. it fires before every sibling
 % has arrived. The sibling count is not re-derived here, so a declaration
 % with k >= n reads as a quorum; use SN_JOIN_QUORUM where the branch count is
 % known and the distinction matters, as the fork-join fixed point does.
 %
 % @par Syntax:
 % @code
 % bool = sn_has_quorum_join(sn)
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
 % <tr><td>bool<td>True if some join declares a positive quorum
 % </table>
%}
function bool = sn_has_quorum_join(sn)

bool = false;
if ~isfield(sn,'nodeparam') || ~iscell(sn.nodeparam)
    return
end
for ind = 1:numel(sn.nodeparam)
    np = sn.nodeparam{ind};
    if ~isstruct(np) || ~isfield(np,'joinStrategy') || ~isfield(np,'joinRequired')
        continue
    end
    for r = 1:numel(np.joinStrategy)
        if isempty(np.joinStrategy{r}) || np.joinStrategy{r} == JoinStrategy.STD
            continue
        end
        if r <= numel(np.joinRequired) && ~isempty(np.joinRequired{r}) && np.joinRequired{r} > 0
            bool = true;
            return
        end
    end
end
end
