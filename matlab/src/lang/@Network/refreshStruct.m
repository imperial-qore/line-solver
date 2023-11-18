function refreshStruct(self, hardRefresh)
% REFRESHSTRUCT()
%
% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

sanitize(self);
if nargin<2
    hardRefresh = true;
end

%% store invariant information
if self.hasStruct && ~hardRefresh
    rtorig = self.sn.rtorig; % this must be destroyed with resetNetwork
end

if self.hasStruct && ~hardRefresh
    nodetypes = sn.nodetypes;
    classnames = sn.classnames;
    nodenames = sn.nodenames;
    refstat = sn.refstat;
else
    nodetypes = getNodeTypes(self);
    classnames = getClassNames(self);
    nodenames = getNodeNames(self);
    refstat = getReferenceStations(self);
end
conn = self.getConnectionMatrix;
njobs = getNumberOfJobs(self);
numservers = getStationServers(self);
lldscaling = getLimitedLoadDependence(self);
cdscaling = getLimitedClassDependence(self);

%% init minimal structure
sn = NetworkStruct(); % create in self to ensure propagation
if isempty(self.sn)
    sn.rtorig = [];
else
    sn.rtorig = self.sn.rtorig;
end
sn.nnodes = numel(nodenames);
sn.nclasses = length(classnames);

%% get routing strategies
routing = zeros(sn.nnodes, sn.nclasses);
for ind=1:sn.nnodes
    for r=1:sn.nclasses
        if isempty(self.nodes{ind}.output.outputStrategy{r})
            routing(ind,r) = RoutingStrategy.ID_DISABLED;
        else
            routing(ind,r) = RoutingStrategy.toId(self.nodes{ind}.output.outputStrategy{r}{2});
        end
    end
end
sn.isslc = false(sn.nclasses,1);
for r=1:sn.nclasses
    if isa(self.classes{r},'SelfLoopingClass')
        sn.isslc(r) = true;
    end
end
sn.nclosedjobs = sum(njobs(isfinite(njobs)));
sn.nservers = numservers;
sn.isstation = NodeType.isStation(nodetypes);
sn.nstations = sum(sn.isstation);
sn.nodetype = -1*ones(sn.nstations,1);
sn.scv = ones(sn.nstations,sn.nclasses);
sn.njobs = njobs(:)';
sn.refstat = refstat;
sn.space = cell(sn.nstations,1);
sn.routing = routing;
sn.chains = [];
sn.lst = {};
sn.lldscaling = lldscaling;
sn.cdscaling = cdscaling;
sn.nodetype = nodetypes;
sn.nstations = sum(sn.isstation);
sn.isstateful = NodeType.isStateful(nodetypes);
sn.isstatedep = false(sn.nnodes,3); % col 1: buffer, col 2: srv, col 3: routing
sn.isfunction = [];
for ind = 1:sn.nstations
    if isa(self.stations{ind},'Queue')
        sn.isfunction(ind) = ~isempty(self.stations{ind}.setupTime);
    end
end
for ind=1:sn.nnodes
    switch sn.nodetype(ind)
        case NodeType.Cache
            sn.isstatedep(ind,2) = true; % state dependent service
            %                            sn.isstatedep(ind,3) = true; % state dependent routing
            %         case NodeType.Place
            %             self.nodes{ind}.init();
            %         case NodeType.Transition
            %             self.nodes{ind}.init(); % this erases enablingConditions
    end

    for r=1:sn.nclasses
        switch sn.routing(ind,r)
            case {RoutingStrategy.ID_RROBIN, RoutingStrategy.ID_WRROBIN, RoutingStrategy.ID_JSQ, RoutingStrategy.ID_RL}
                sn.isstatedep(ind,3) = true; % state dependent routing
        end
    end
end
sn.nstateful = sum(sn.isstateful);
sn.state = cell(sn.nstations,1);
for i=1:sn.nstateful
    sn.state{i} = [];
end
sn.nodenames = nodenames;
sn.classnames = classnames;
sn.connmatrix = conn;
sn.nodeparam{ind}.modenames = cell(sn.nnodes,1);
sn.nodeparam{ind}.nmodes = zeros(sn.nnodes,1);
sn.nodeparam{ind}.nmodeservers = cell(sn.nnodes,1);
sn.nodeparam{ind}.enabling = cell(sn.nnodes,1);
sn.nodeparam{ind}.firing = cell(sn.nnodes,1);
sn.nodeparam{ind}.firingprio = cell(sn.nnodes,1);
sn.nodeparam{ind}.fireweight = cell(sn.nnodes,1);
sn.nodeparam{ind}.firingid = cell(sn.nnodes,1);
sn.nodeparam{ind}.firingproc = cell(sn.nnodes,1);
sn.nodeparam{ind}.firingprocid = cell(sn.nnodes,1);
sn.nodeparam{ind}.firingphases = cell(sn.nnodes,1);
sn.nodeparam{ind}.inhibiting = cell(sn.nnodes,1);

sn.nodeToStateful =[];
sn.nodeToStation =[];
sn.stationToNode =[];
sn.stationToStateful =[];
sn.statefulToNode =[];
for ind=1:sn.nnodes
    sn.nodeToStateful(ind) = nd2sf(sn,ind);
    sn.nodeToStation(ind) = nd2st(sn,ind);
end
for ist=1:sn.nstations
    sn.stationToNode(ist) = st2nd(sn,ist);
    sn.stationToStateful(ist) = st2sf(sn,ist);
end
for isf=1:sn.nstateful
    sn.statefulToNode(isf) = sf2nd(sn,isf);
end
sn.fj = self.getForkJoins();
self.sn = sn;
refreshPriorities(self);
refreshService(self);
if any(nodetypes == NodeType.Cache)
    refreshChains(self, false); % wantVisits
else
    refreshChains(self, true); % wantVisits
end
sn = self.sn;
refclasses = getReferenceClasses(self);
refclass = zeros(1,sn.nchains);
for c=1:sn.nchains
    isect = intersect(sn.inchain{c},find(refclasses));
    if any(isect)
        refclass(c) = isect;
    end
end
sn.refclass = refclass;
self.sn = sn;

refreshLocalVars(self); % depends on chains (rtnodes)
refreshSync(self); % this assumes that refreshChain is called before
refreshPetriNetNodes(self);
self.hasStruct = true;
end

function stat_idx = nd2st(sn, node_idx)
% STAT_IDX = ND2ST(NODE_IDX)

if sn.isstation(node_idx)
    stat_idx = at(cumsum(sn.isstation),node_idx);
else
    stat_idx = NaN;
end
end

function node_idx = st2nd(sn,stat_idx)
% NODE_IDX = ST2ND(SELF,STAT_IDX)

v = cumsum(sn.isstation) == stat_idx;
if any(v)
    node_idx =  find(v, 1);
else
    node_idx = NaN;
end
end

function sful_idx = st2sf(sn,stat_idx)
% SFUL_IDX = ST2SF(SELF,STAT_IDX)

sful_idx = nd2sf(sn,st2nd(sn,stat_idx));
end

function sful_idx = nd2sf(sn, node_idx)
% SFUL_IDX = ND2SF(NODE_IDX)

if sn.isstateful(node_idx)
    sful_idx = at(cumsum(sn.isstateful),node_idx);
else
    sful_idx = NaN;
end
end

function node_idx = sf2nd(sn,stat_idx)
% NODE_IDX = SF2ND(SELF,STAT_IDX)

v = cumsum(sn.isstateful) == stat_idx;
if any(v)
    node_idx =  find(v, 1);
else
    node_idx = NaN;
end
end
