function plotGraphSimple(self, method)
% PLOTGRAPHSIMPLE(SELF, METHOD)
%
% Plot graph without colors, suitable for inclusion in scientific papers
% METHOD: hashnames, names, hashids, ids

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<2 %~exist('useNodes','var')
    method = 'nodes';
end

lqn = self.getStruct;
T = lqn.graph;
figure;
switch method
    case {'nodes','hashnames'}
        lqn.hashnames = strrep(lqn.hashnames,'\_','_'); % remove escaping if it exists
        lqn.hashnames = strrep(lqn.hashnames,'_','\_'); % reapply it include to those which did not have it
        h = plot(digraph(T),'Layout','layered','NodeLabel',lqn.hashnames);
    case 'names'
        lqn.names = strrep(lqn.names,'\_','_'); % remove escaping if it exists
        lqn.names = strrep(lqn.names,'_','\_'); % reapply it include to those which did not have it
        h = plot(digraph(T),'Layout','layered','NodeLabel',lqn.names);
    case 'ids'
        h = plot(digraph(T),'Layout','layered');
    case 'hashids'
        hashids = lqn.hashnames;
        for i=1:length(lqn.hashnames)
            hashids{i} = [lqn.hashnames{i}(1:2),num2str(i)];
        end
        h = plot(digraph(T),'Layout','layered','NodeLabel',hashids);
end
title(['Model: ',self.name]);

%for r=find(lqnGraph.Edges.Pre==1)' %AndFork
%    highlight(h,lqnGraph.Edges.EndNodes(r,:),'EdgeColor','m')
%end
%for r=find(lqnGraph.Edges.Post==1)' % AndJoin
%    highlight(h,lqnGraph.Edges.EndNodes(r,:),'EdgeColor','m')
%end

for r=find(lqn.type==LayeredNetworkElement.HOST)'
    if r>0
        highlight(h,r,'NodeColor','black');
    end
end

for r=find(lqn.type==LayeredNetworkElement.ACTIVITY)'
    if r>0
        highlight(h,r,'NodeColor','black');
    end
end

for r=find(lqn.type==LayeredNetworkElement.ACTIVITY)'
    p = find(lqn.graph(:,r))';
    if r>0
        switch lqn.actpretype(r)
            case ActivityPrecedenceType.PRE_AND
                % highlight(h,p,'Marker','o');
        end
    end
    p = find(lqn.graph(r,:))';
    if r>0
        switch lqn.actposttype(r)
            case ActivityPrecedenceType.POST_AND
                % highlight(h,p,'Marker','o');
        end
    end
end

for r=find(lqn.type==LayeredNetworkElement.TASK)'
    if r>0
        if lqn.isref(r)
            highlight(h,r,'NodeColor','black');
        else
            highlight(h,r,'NodeColor','black');
        end
    end
end

for r=find(lqn.type==LayeredNetworkElement.ENTRY)'
    if r>0
        highlight(h,r,'NodeColor','black');
    end
end

mult = nan(lqn.nidx,1);
mult(1:(lqn.tshift+lqn.ntasks),1)=lqn.mult;
row = dataTipTextRow('multiplicity',mult);
h.DataTipTemplate.DataTipRows(2) = row;
D = 0*lqn.mult;
for i=1:length(lqn.hostdem)
    if ~isempty(lqn.hostdem{i})
        D(i) = lqn.hostdem{i}.getMean;
    end
end
row = dataTipTextRow('hostDemand',D);
h.DataTipTemplate.DataTipRows(end+1) = row;
sched = cell(lqn.nidx,1);
for i=1:length(sched)
    sched{i} = 'n/a';
end
sched(1:(lqn.tshift+lqn.ntasks),1)=lqn.sched;
row = dataTipTextRow('scheduling',sched);
h.DataTipTemplate.DataTipRows(end+1) = row;
end
