function V=pathsbetween(dag,node1,node2,Vin)
% determines all paths between two nodes on a dag
if nargin < 4
    Vin = [];
end
if node1 == node2
    V = Vin;
    return
end
V = [];
if isempty(find(dag(node1,:)))
    return
end
for successor=find(dag(node1,:))
    Vin = [Vin,successor];
    Vnew = pathsbetween(dag,successor,node2,Vin);
    if ~isempty(Vnew)
        if isempty(V)
            V= Vnew;
        else
            V(end+1,1:length(Vnew)) = Vnew; %[successor,Vnew];
        end
    end
end
end