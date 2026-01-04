function [H,G] = getGraph(self)
% [H,G] = GETGRAPH()

G = digraph(); TG = Table();
M = self.getNumberOfNodes;
K = self.getNumberOfClasses;
sn = self.getStruct;
[P,Pnodes] = getRoutingMatrix(self);
name = {}; sched = {}; type = {}; nservers = [];
for ist=1:M
    name{end+1} = self.nodes{ist}.name;
    type{end+1} = class(self.nodes{ist});
    if ~isa(self.nodes{ist},'Join')
        sched{end+1} = self.nodes{ist}.schedStrategy;
    else
        sched{end+1} = '';
    end
    if isa(self.nodes{ist},'Station')
        nservers(end+1) = self.nodes{ist}.getNumberOfServers;
    else
        nservers(end+1) = 0;
    end
end
TG.Name = name(:);
TG.Type = type(:);
TG.Sched = sched(:);
TG.Servers = nservers(:);
G = G.addnode(TG);
for ist=1:M
    for jst=1:M
        for k=1:K
            if Pnodes((ist-1)*K+k,(jst-1)*K+k) > 0
                G = G.addedge(self.nodes{ist}.name,self.nodes{jst}.name, Pnodes((ist-1)*K+k,(jst-1)*K+k));
            end
        end
    end
end
H = digraph(); TH = Table();
I = self.getNumberOfStations;
name = {}; sched = {}; type = {}; jobs = zeros(I,1); nservers = [];
for ind=1:I
    name{end+1} = self.stations{ind}.name;
    type{end+1} = class(self.stations{ind});
    if ~isa(self.stations{ind},'Join')
        sched{end+1} = self.stations{ind}.schedStrategy;
    else
        sched{end+1} = '';
    end
    for k=1:K
        if sn.refstat(k)==ind
            jobs(ind) = jobs(ind) + sn.njobs(k);
        end
    end
    if isa(self.nodes{ind},'Station')
        nservers(end+1) = self.nodes{ind}.getNumberOfServers;
    else
        nservers(end+1) = 0;
    end
end
TH.Name = name(:);
TH.Type = type(:);
TH.Sched = sched(:);
TH.Jobs = jobs(:);
TH.Servers = nservers(:);
H = H.addnode(TH);
rate = [];
classes = {};
for ind=1:I
    for jnd=1:I
        for k=1:K
            if P((ind-1)*K+k,(jnd-1)*K+k) > 0
                rate(end+1) = sn.rates(ind,k);
                classes{end+1} = self.classes{k}.name;
                H = H.addedge(self.stations{ind}.name, self.stations{jnd}.name, P((ind-1)*K+k,(jnd-1)*K+k));
            end
        end
    end
end
H.Edges.Rate = rate(:);
H.Edges.Class = classes(:);
H = H.rmedge(find(isnan(H.Edges.Rate)));
sourceObj = self.getSource;
if ~isempty(sourceObj)
    %                 sink = self.getSink;
    %                 H=H.addnode(sink.name);
    %                 H.Nodes.Type{end}='Sink';
    %                 H.Nodes.Sched{end}='ext';
    %H = H.rmedge(find(isnan(H.Edges.Rate)));
    %sourceIdx = model.getIndexSourceNode;
    %                toDel = findstring(H.Edges.EndNodes(:,2),sourceObj.name);
    %                for j=toDel(:)'
    %                    H = H.rmedge(j);
    %                end
end
end