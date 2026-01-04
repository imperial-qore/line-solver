%{ @file lsn_max_multiplicity.m
 %  @brief Computes the maximum multiplicity of nodes in a Layered Software Network
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes maximum multiplicity (concurrency level) for LSN nodes
 %
 % @details
 % This function computes the maximum multiplicity (throughput capacity) for
 % each node in a Layered Software Network (LSN). Uses Kahn's algorithm for
 % topological sorting to propagate constraints through the network.
 %
 % @par Syntax:
 % @code
 % outflow = lsn_max_multiplicity(lsn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lsn<td>Layered Software Network structure
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>outflow<td>Maximum multiplicity (throughput capacity) for each node
 % </table>
%}
function outflow = lsn_max_multiplicity(lsn)
ag = lsn.dag > 0;
mult = lsn.mult;
type = lsn.type;
isref = lsn.isref;
n = size(ag, 1);

% Manual topological sort (Kahn's algorithm)
order = kahn(ag);

% initially load ref task multiplicity
inflow = zeros(n,1);
for i = 1:n
    if type(i) == LayeredNetworkElement.TASK && isref(i)
        inflow(i) = mult(i);
    % Also account for entries with open arrivals
    elseif type(i) == LayeredNetworkElement.ENTRY
        if isfield(lsn, 'arrival') && ~isempty(lsn.arrival) && ...
           iscell(lsn.arrival) && i <= length(lsn.arrival) && ...
           ~isempty(lsn.arrival{i})
            % Entry has open arrival - needs at least 1 thread of its parent task
            inflow(i) = 1;
        end
    end
end
outflow = zeros(n,1);

if length(mult) < n
    mult(end+1:n) = Inf;
end

for k = 1:n
    i = order(k);
    outflow(i) = min(inflow(i), mult(i));
    for j = [1:i-1,i+1:n]
        if ag(i,j)
            inflow(j) = inflow(j) + outflow(i);
        end
    end
end
%inflow(type > LayeredNetworkElement.TASK )=0;
for i = 1:n
    if type(i) == LayeredNetworkElement.TASK && mult(i)==Inf && ~isref(i)
        outflow(i) = Inf;
    end
end
end
