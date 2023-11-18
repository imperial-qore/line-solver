function self = linkFromNodeRoutingMatrix(self, P)
% SELF = LINK(P)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

sanitize(self);

isReset = false;
if ~isempty(self.sn)
    isReset = true;
    self.resetNetwork; % remove artificial class switch nodes
end
R = self.getNumberOfClasses;
Mnodes = self.getNumberOfNodes;

if iscell(P) && R>1
    line_error(mfilename,'Multiclass model: the linked routing matrix P must be a matrix');
end

for r=1:R
    for s=1:R
        Pnodes{r,s} = zeros(R);
        for i=1:Mnodes
            for j=1:Mnodes
                if P((i-1)*R+r,(j-1)*R+s)>0
                    Pnodes{r,s}(i,j) = P((i-1)*R+r,(j-1)*R+s);
                end
            end
        end
    end
end

connected = zeros(Mnodes);
nodes = self.nodes;
for r=1:R
    [I,J,S] = find(Pnodes{r,r});
    for k=1:length(I)        
        if connected(I(k),J(k)) == 0
            self.addLink(nodes{I(k)}, nodes{J(k)});
            connected(I(k),J(k)) = 1;
        end
        nodes{I(k)}.setProbRouting(self.classes{r}, nodes{J(k)}, S(k));
    end
end
self.nodes = nodes;

for i=1:Mnodes
    if isa(self.nodes{i},'Place')
        self.nodes{i}.init;
    end
end

if isReset
    self.refreshChains; % without this exception with linkAndLog
end

end
