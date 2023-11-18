function [H,G] = getGraph(self)
% [H,G] = GETGRAPH()

G = digraph(); TG = Table();
M = self.getNumberOfNodes;
K = self.getNumberOfClasses;
sn = self.getStruct;
[P,Pnodes] = getRoutingMatrix(self);
name = {}; sched = {}; type = {}; nservers = [];
for i=1:M
    name{end+1} = self.nodes{i}.name;
    type{end+1} = class(self.nodes{i});
    if ~isa(self.nodes{i},'Join')
        sched{end+1} = self.nodes{i}.schedStrategy;
    else
        sched{end+1} = '';
    end
    if isa(self.nodes{i},'Station')
        nservers(end+1) = self.nodes{i}.getNumberOfServers;
    else
        nservers(end+1) = 0;
    end
end
TG.Name = name(:);
TG.Type = type(:);
TG.Sched = sched(:);
TG.Servers = nservers(:);
G = G.addnode(TG);
for i=1:M
    for j=1:M
        for k=1:K
            if Pnodes((i-1)*K+k,(j-1)*K+k) > 0
                G = G.addedge(self.nodes{i}.name,self.nodes{j}.name, Pnodes((i-1)*K+k,(j-1)*K+k));
            end
        end
    end
end
H = digraph(); TH = Table();
I = self.getNumberOfStations;
name = {}; sched = {}; type = {}; jobs = zeros(I,1); nservers = [];
for i=1:I
    name{end+1} = self.stations{i}.name;
    type{end+1} = class(self.stations{i});
    if ~isa(self.stations{i},'Join')
        sched{end+1} = self.stations{i}.schedStrategy;
    else
        sched{end+1} = '';
    end
    for k=1:K
        if sn.refstat(k)==i
            jobs(i) = jobs(i) + sn.njobs(k);
        end
    end
    if isa(self.nodes{i},'Station')
        nservers(end+1) = self.nodes{i}.getNumberOfServers;
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
for i=1:I
    for j=1:I
        for k=1:K
            if P((i-1)*K+k,(j-1)*K+k) > 0
                rate(end+1) = sn.rates(i,k);
                classes{end+1} = self.classes{k}.name;
                H = H.addedge(self.stations{i}.name, self.stations{j}.name, P((i-1)*K+k,(j-1)*K+k));
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