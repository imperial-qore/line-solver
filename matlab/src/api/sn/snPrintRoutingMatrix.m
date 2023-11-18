function snPrintRoutingMatrix(sn, onlyclass)
% SNPRINTROUTINGMATRIX()

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
                    if sn.nodetype(i) == NodeType.ID_CACHE
                        pr = 'state-dependent';
                    elseif sn.nodetype(i) == NodeType.ID_SINK
                        continue
                    else 
                        if sn.routing(i,r) == RoutingStrategy.ID_DISABLED
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
end