%{ @file sn_rtnodes_to_rtorig.m
 %  @brief Converts node routing matrix to the original routing matrix format
 %
 %  @author LINE Development Team
%}

%{
 % @brief Converts node routing matrix to the original routing matrix format
 %
 % @details
 % This function converts the node-level routing matrix to the original
 % routing matrix format, excluding class-switching nodes.
 %
 % @par Syntax:
 % @code
 % [rtorigcell, rtorig] = sn_rtnodes_to_rtorig(sn)
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
 % <tr><td>rtorigcell<td>Cell array representation of the routing matrix
 % <tr><td>rtorig<td>Sparse matrix representation of the routing matrix
 % </table>
%}
function [rtorigcell,rtorig] = sn_rtnodes_to_rtorig(sn)
K = sn.nclasses;
rtnodes = sn.rtnodes;

csshift = sn.nnodes;
for ind = 1:sn.nnodes
    if startsWith(sn.nodenames{ind}, 'CS_')
        csshift = ind-1;
        break
    end
end

colToKeep=[];
for ind = 1:csshift
    for k = 1:K
        colToKeep(end+1) = (ind-1)*K+k;
    end
end

rtorig = dtmc_stochcomp(rtnodes, colToKeep);
rtorigcell = cellzeros(K,K,csshift,csshift);

% for cache rt, replace NaNs for unknown probabilities with 0
rtorig(isnan(rtorig)) = 0;

for ind = 1:csshift
    if sn.nodetype(ind) ~= NodeType.Sink
        for jnd = 1:csshift
            for r = 1:K
                for s = 1:K
                    rtorigcell{r,s}(ind,jnd) = rtorig((ind-1)*K+r,(jnd-1)*K+s);
                end
            end
        end
    end
end

end