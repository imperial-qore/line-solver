%{ @file sn_print_routing_matrix.m
 %  @brief Prints the routing matrix of the network
 %
 %  @author LINE Development Team
%}

%{
 % @brief Prints the routing matrix of the network
 %
 % @details
 % This function displays the routing probabilities between nodes and classes
 % in a human-readable format.
 %
 % @par Syntax:
 % @code
 % sn_print_routing_matrix(sn)
 % sn_print_routing_matrix(sn, onlyclass)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>onlyclass<td>(Optional) Filter printing for a specific class
 % </table>
%}
function sn_print_routing_matrix(sn, onlyclass)

node_names = sn.nodenames;
classnames = sn.classnames;
rtnodes = sn.rtnodes;
nnodes = sn.nnodes;
nclasses = sn.nclasses;
for i=1:nnodes
    for r=1:nclasses
        for j=1:nnodes
            for s=1:nclasses
                if rtnodes((i-1)*nclasses+r,(j-1)*nclasses+s)>0
                    if sn.nodetype(i) == NodeType.Cache
                        pr = 'state-dependent';
                    elseif sn.nodetype(i) == NodeType.Sink
                        continue
                    else 
                        if sn.routing(i,r) == RoutingStrategy.DISABLED
                            %pr = 'Disabled';
                            continue
                        else
                            pr = num2str(rtnodes((i-1)*nclasses+r,(j-1)*nclasses+s),'%f');
                        end
                    end
                    if nargin==1
                        line_printf('\n%s [%s] => %s [%s] : Pr=%s',node_names{i}, classnames{r}, node_names{j}, classnames{s}, pr);
                    else
                        if strcmpi(classnames{r},onlyclass.name) || strcmpi(classnames{s},onlyclass.name)
                            line_printf('\n%s [%s] => %s [%s] : Pr=%s',node_names{i}, classnames{r}, node_names{j}, classnames{s}, pr);
                        end
                    end
                end
            end
        end
    end
end
line_printf('\n');
end