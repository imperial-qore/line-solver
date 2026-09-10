function self = link(self, P)
% SELF = LINK(P)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isempty(self.connections)
    line_error(mfilename,'The Network.link method cannot be used after calling the addLink() method. Use Node.setProbRouting instead to configure routing probabilities.');
end

if isa(P,'RoutingMatrix')
    P = P.getCell();
end
sanitize(self);

isReset = false;
if ~isempty(self.sn)
    isReset = true;
    self.resetNetwork; % remove artificial class switch nodes
end
K = self.getNumberOfClasses;
I = self.getNumberOfNodes;

if ~iscell(P) && K>1
    line_error(mfilename,'Multiclass model: the linked routing matrix P must be a cell array, e.g., P = model.initRoutingMatrix; P{1} = Pclass1; P{2} = Pclass2.');
end

isLinearP = true;
if size(P,1) == size(P,2)
    for s=2:K
        for r=1:K
            if nnz(P{r,s})>0
                isLinearP = false;
                break;
            end
        end
    end
    % P may look linear only because state-dependent routing leaves some zero entries unspecified
    %cacheNodes = find(cellfun(@(c) isa(c,'Cache'), self.getStatefulNodes));
    for ind=1:I
        switch class(self.nodes{ind})
            case 'Cache'
                % A cache needs class-switch to distinguish hits/misses unless the model is degenerate
                isLinearP = false;
                % Both are sparse and indexed by class, and a retrieval system extends
                % missClass past hitClass, so compare them zero-padded to one length
                hitCls_ = full(self.nodes{ind}.server.hitClass(:)).';
                missCls_ = full(self.nodes{ind}.server.missClass(:)).';
                L_ = max(length(hitCls_), length(missCls_));
                hitCls_(end+1:L_) = 0;
                missCls_(end+1:L_) = 0;
                if L_ > 0 && all(hitCls_ == missCls_)
                    line_warning(mfilename,'Ambiguous use of hitClass and missClass at cache, it is recommended to use different classes.\n');
                end
        end
    end
end

% This block is to make sure that P = model.initRoutingMatrix; P{2} writes
% into P{2,2} rather than being interpreted as P{2,1}.
if isLinearP
    Ptmp = P;
    P = cell(K,K);
    for r=1:K
        if iscell(Ptmp)
            P{r,r} = Ptmp{r};
        else
            P{r,r} = Ptmp;
        end
        for s=1:K
            if s~=r
                P{r,s} = 0*Ptmp{r};
            end
        end
    end
end

% assign routing for self-looping jobs
for r=1:K
    if isa(self.classes{r},'SelfLoopingClass')
        for s=1:K
            P{r,s} = 0 * P{r,s};
        end
        P{r,r}(self.classes{r}.refstat, self.classes{r}.refstat) = 1.0;
    end
end

% link virtual sinks automatically to sink
ispool = cellisa(self.nodes,'Sink');
if sum(ispool) > 1
    line_error(mfilename,'The model can have at most one sink node.');
end

if sum(cellisa(self.nodes,'Source')) > 1
    line_error(mfilename,'The model can have at most one source node.');
end
ispool_nnz = find(ispool)';


if ~iscell(P)
    if K>1
        newP = cell(1,K);
        for r=1:K
            newP{r} = P;
        end
        P = newP;
    else %R==1
        % single class
        for ind=ispool_nnz
            P((ind-1)*K+1:ind*K,:)=0;
        end
        Pmat = P;
        P = cell(K,K);
        for r=1:K
            for s=1:K
                P{r,s} = zeros(I);
                for ind=1:I
                    for jnd=1:I
                        P{r,s}(ind,jnd) = Pmat((ind-1)*K+r,(jnd-1)*K+s);
                    end
                end
            end
        end
    end
end

if numel(P) == K
    % 1 matrix per class
    for r=1:K
        for ind=ispool_nnz
            P{r}((ind-1)*K+1:ind*K,:)=0;
        end
    end
    Pmat = P;
    P = cell(K,K);
    for r=1:K
        P{r,r} = Pmat{r};
        for s=setdiff(1:K,r)
            P{r,s} = zeros(I);
        end
    end
end

% Default per-item retrieval routing inherited from the read class -- see _kb/09-ldes-and-cache.md
for ind=1:I
    if isa(self.nodes{ind}, 'Cache') && isprop(self.nodes{ind}, 'retrievalSystemQueueIndices') ...
            && ~isempty(keys(self.nodes{ind}.retrievalSystemQueueIndices))
        cacheNode = self.nodes{ind};
        c = cacheNode.index;
        nItems = cacheNode.items.nitems;
        rcMap = cacheNode.retrievalSystemQueueIndices;
        ks = keys(rcMap);
        for kk=1:numel(ks)
            r = double(ks(kk)) + 1;      % read class index (1-based; key is index-1)
            Q = rcMap{ks(kk)};           % retrieval-system queue node indices
            nQ = numel(Q);
            if isempty(P{r,r})
                continue
            end
            Pr = P{r,r};
            for it=1:nItems
                rClass = cacheNode.server.retrievalClasses(it, r);
                if rClass <= 0
                    continue
                end
                if isempty(P{rClass,rClass})
                    P{rClass,rClass} = zeros(I);
                end
                for q=1:nQ
                    if Pr(c, Q(q)) > 0                 % cache -> queue (entry)
                        P{rClass,rClass}(c, Q(q)) = Pr(c, Q(q));
                    end
                    if Pr(Q(q), c) > 0                 % queue -> cache (exit)
                        P{rClass,rClass}(Q(q), c) = Pr(Q(q), c);
                    end
                    for dq=1:nQ                        % queue -> queue
                        if Pr(Q(q), Q(dq)) > 0
                            P{rClass,rClass}(Q(q), Q(dq)) = Pr(Q(q), Q(dq));
                        end
                    end
                end
            end
            % consume the read class's template edges over the queue set
            for q=1:nQ
                P{r,r}(c, Q(q)) = 0;
                P{r,r}(Q(q), c) = 0;
                for dq=1:nQ
                    P{r,r}(Q(q), Q(dq)) = 0;
                end
            end
        end
    end
end

% Inject deferred retrieval-system routing entries from Cache.setRetrievalSystem -- see _kb/09-ldes-and-cache.md
for ind=1:I
    if isa(self.nodes{ind}, 'Cache') && isprop(self.nodes{ind}, 'retrievalRoutingEntries') ...
            && ~isempty(self.nodes{ind}.retrievalRoutingEntries)
        rre = self.nodes{ind}.retrievalRoutingEntries;
        for e=1:numel(rre)
            ent = rre{e}; % [fromCls, toCls, srcNode, dstNode, prob]
            if isempty(P{ent(1),ent(2)})
                P{ent(1),ent(2)} = zeros(I);
            end
            P{ent(1),ent(2)}(ent(3),ent(4)) = ent(5);
        end
    end
end

isemptyP = false(K,K);
for r=1:K
    for s=1:K
        if isempty(P{r,s})
            isemptyP(r,s)= true;
            P{r,s} = zeros(I);
        else
            for ind=ispool_nnz
                P{r,s}(ind,:)=0;
            end
        end
    end
end

csnodematrix = cell(I,I);
for ind=1:I
    for jnd=1:I
        csnodematrix{ind,jnd} = zeros(K,K);
    end
end

for r=1:K
    for s=1:K
        if ~isemptyP(r,s)
            [If,Jf] = find(P{r,s});
            for k=1:size(If,1)
                csnodematrix{If(k),Jf(k)}(r,s) = P{r,s}(If(k),Jf(k));
            end
        end
    end
end


%             for r=1:R
%                 Psum=cellsum({P{r,:}})*ones(M,1);
%                 if min(Psum)<1-GlobalConstants.CoarseTol
%                   line_error(mfilename,'Invalid routing probabilities (Node %d departures, switching from class %d).',minpos(Psum),r);
%                 end
%                 if max(Psum)>1+GlobalConstants.CoarseTol
%                   line_error(mfilename,sprintf('Invalid routing probabilities (Node %d departures, switching from class %d).',maxpos(Psum),r));
%                 end
%             end

self.sn.rtorig = P;

% As we will now create a CS for each link i->j,
% we now condition on the job going from node i to j
for ind=1:I
    for jnd=1:I
        for r=1:K
            S = sum(csnodematrix{ind,jnd}(r,:));
            if S>0
                csnodematrix{ind,jnd}(r,:)=csnodematrix{ind,jnd}(r,:)/S;
            else
                csnodematrix{ind,jnd}(r,r)=1.0;
            end
        end
    end
end

csid = zeros(I);
% eye(K): a job to an autoAdded class switch stays in the same class beforehand; csMatrix diagonal is irrelevant for chains
csMatrix = eye(K);
nodeNames = self.getNodeNames;
for ind=1:I
    for jnd=1:I
        csMatrix = csMatrix + csnodematrix{ind,jnd};
        if ~isdiag(csnodematrix{ind,jnd})
            self.nodes{end+1} = ClassSwitch(self, sprintf('CS_%s_to_%s',nodeNames{ind},nodeNames{jnd}),csnodematrix{ind,jnd});
            self.nodes{end}.autoAdded = true;
            csid(ind,jnd) = length(self.nodes);
        end
    end
end

for ind=1:I
    % this is to ensure that also stateful cs like caches
    % are accounted
    if isa(self.nodes{ind},'Cache')
        for r=find(self.nodes{ind}.server.hitClass)
            csMatrix(r,self.nodes{ind}.server.hitClass(r)) = 1.0;
        end
        for r=find(self.nodes{ind}.server.missClass)
            csMatrix(r,self.nodes{ind}.server.missClass(r)) = 1.0;
        end
    elseif isa(self.nodes{ind},'ClassSwitch')
        if isempty(self.nodes{ind}.server.csMatrix )
            line_error(mfilename,'Uninitialized ClassSwitch node, use the setClassSwitchingMatrix method.');
        end
        csMatrix  = csMatrix | self.nodes{ind}.server.csMatrix > 0.0;
    end
end

self.csMatrix = csMatrix~=0;

Ip = length(self.nodes); % number of nodes after addition of cs nodes

% resize matrices
for r=1:K
    for s=1:K
        P{r,s}((I+1):Ip,(I+1):Ip)=0;
    end
end

for ind=1:I
    for jnd=1:I
        if csid(ind,jnd)>0
            % re-route
            for r=1:K
                for s=1:K
                    if P{r,s}(ind,jnd)>0
                        P{r,r}(ind,csid(ind,jnd)) = P{r,r}(ind,csid(ind,jnd)) + P{r,s}(ind,jnd);
                        P{r,s}(ind,jnd) = 0;
                        P{s,s}(csid(ind,jnd),jnd) = 1;
                    end
                end
            end
        end
    end
end

connected = zeros(Ip);
nodes = self.nodes;
% Clear non-station nodes' outputStrategy before setting new routing (prevents accumulation across link() calls) -- see _kb/04-networkstruct.md
for ind=1:Ip
    if ~isa(nodes{ind}, 'Station') && ismethod(nodes{ind}.output, 'initDispatcherJobClasses')
        if ismethod(nodes{ind}.output, 'initDispatcherJobClassesPreservingRouted')
            nodes{ind}.output.initDispatcherJobClassesPreservingRouted(self.classes);
        else
            nodes{ind}.output.initDispatcherJobClasses(self.classes);
        end
        % Sync sn.routing with the post-clear outputStrategy so getRoutingMatrix won't try PROB routing for cleared classes
        if ~isempty(self.sn) && isfield(self.sn, 'routing') && ind <= size(self.sn.routing, 1)
            for k=1:K
                if k <= numel(nodes{ind}.output.outputStrategy) && ~isempty(nodes{ind}.output.outputStrategy{k})
                    entry = nodes{ind}.output.outputStrategy{k};
                    self.sn.routing(ind,k) = RoutingStrategy.fromText(entry{2});
                else
                    self.sn.routing(ind,k) = RoutingStrategy.DISABLED;
                end
            end
        end
    end
end
for r=1:K
    [If,Jf,S] = find(P{r,r});
    for k=1:length(If)
        if connected(If(k),Jf(k)) == 0
            self.addLink(nodes{If(k)}, nodes{Jf(k)});
            connected(If(k),Jf(k)) = 1;
        end
        nodes{If(k)}.setProbRouting(self.classes{r}, nodes{Jf(k)}, S(k));
    end
end
self.nodes = nodes;

% Refresh sn.routing from outputStrategy: initDispatcherJobClasses above set it DISABLED, setProbRouting has since repopulated it
if ~isempty(self.sn) && isfield(self.sn, 'routing')
    for ind=1:Ip
        if ~isa(nodes{ind}, 'Station') && ind <= size(self.sn.routing, 1)
            for k=1:K
                if isempty(nodes{ind}.output.outputStrategy{k})
                    self.sn.routing(ind,k) = RoutingStrategy.DISABLED;
                else
                    self.sn.routing(ind,k) = RoutingStrategy.fromText(nodes{ind}.output.outputStrategy{k}{2});
                end
            end
        end
    end
end

% Validate the routing probabilities -- see _kb/04-networkstruct.md
% Pcheck normalizes the two accepted shapes (plain matrix, or K x K class-pair cell) so both are checked the same way.
if iscell(P)
    Pcheck = P;
else
    Pcheck = {P};
end

% (1) Nonnegativity, checked separately from the row-sum test below -- see _kb/04-networkstruct.md
if self.enableChecks
    for r=1:size(Pcheck,1)
        for s=1:size(Pcheck,2)
            if ~isempty(Pcheck{r,s})
                [Ineg,Jneg] = find(Pcheck{r,s} < -GlobalConstants.FineTol);
                if ~isempty(Ineg)
                    if size(Pcheck,1) > 1 || size(Pcheck,2) > 1
                        line_error(mfilename,sprintf('Negative routing probability %g from node %s to node %s (class %s to class %s). Routing probabilities must be nonnegative.', full(Pcheck{r,s}(Ineg(1),Jneg(1))), self.nodes{Ineg(1)}.name, self.nodes{Jneg(1)}.name, self.classes{r}.name, self.classes{s}.name));
                    else
                        line_error(mfilename,sprintf('Negative routing probability %g from node %s to node %s. Routing probabilities must be nonnegative.', full(Pcheck{r,s}(Ineg(1),Jneg(1))), self.nodes{Ineg(1)}.name, self.nodes{Jneg(1)}.name));
                    end
                end
            end
        end
    end
end

% (2) Row-sum <= 1, summed over destination j AND arrival class s (class-switching correctness) -- see _kb/04-networkstruct.md
if self.enableChecks
    for ind=1:I
        % Place/Transition (SPN incidence arcs) and Router (connectivity declared before setRouting) are exempt -- see _kb/04-networkstruct.md
        if isa(self.nodes{ind},'Place') || isa(self.nodes{ind},'Transition') || isa(self.nodes{ind},'Router')
            continue
        end
        % FORK legitimately emits on several output links at once; guarded since not every node carries schedStrategy (e.g. Join)
        if isprop(self.nodes{ind},'schedStrategy') && ~isempty(self.nodes{ind}.schedStrategy) && SchedStrategy.toId(self.nodes{ind}.schedStrategy) == SchedStrategy.FORK
            continue
        end
        for r=1:size(Pcheck,1)
            pOut = 0;
            for s=1:size(Pcheck,2)
                if ~isempty(Pcheck{r,s})
                    pOut = pOut + sum(Pcheck{r,s}(ind,:));
                end
            end
            if pOut > 1.0 + GlobalConstants.FineTol
                if size(Pcheck,1) > 1
                    line_error(mfilename,sprintf('The total routing probability for jobs leaving node %s in class %s is %g, which is greater than 1.0.',self.nodes{ind}.name,self.classes{r}.name,full(pOut)));
                else
                    line_error(mfilename,sprintf('The total routing probability for jobs leaving node %s is %g, which is greater than 1.0.',self.nodes{ind}.name,full(pOut)));
                end
            end
        end
        %        elseif pOut < 1.0 - GlobalConstants.FineTol % we cannot check this case as class r may not reach station i, in which case its outgoing routing prob is zero
        %            if self.nodes{i}.schedStrategy ~= SchedStrategy.EXT % if not a sink
        %                line_error(mfilename,'The total routing probability for jobs leaving node %s in class %s is less than 1.0.',self.nodes{i}.name,self.classes{r}.name);
        %            end
    end
end

for ind=1:I
    if isa(self.nodes{ind},'Place')
        self.nodes{ind}.init;
    end
end

if isReset && ~isempty(self.sn) && isfield(self.sn,'rates')
    self.refreshChains; % without this exception with linkAndLog
end

%% Check for reducible routing (absorbing states)
if self.enableChecks
    % Pass routing matrix directly to avoid getStruct() call during link()
    [isErg, ergInfo] = self.isRoutingErgodic(self.sn.rtorig);
    if ~isErg && ~isempty(ergInfo.absorbingStations)
        % Build warning message
        absNames = strjoin(ergInfo.absorbingStations, ', ');
        line_warning(mfilename, 'Reducible network topology detected, results may be unreliable.\n');
    end
end

%% Check that order-independent (OI) stations have a permutation-invariant rate
hasOI = false;
for ind = 1:I
    nd = self.nodes{ind};
    if isa(nd, 'Queue') && SchedStrategy.toId(nd.schedStrategy) == SchedStrategy.OI ...
            && ~isempty(nd.svcRateFun)
        hasOI = true; break
    end
end
if self.enableChecks && hasOI
    K = self.getNumberOfClasses;
    pop = zeros(1, K);
    for r = 1:K
        if isa(self.classes{r}, 'ClosedClass')
            pop(r) = self.classes{r}.population;
        else
            pop(r) = Inf;   % open classes have infinite population
        end
    end
    % The conserved quantity under class switching is the CHAIN population, not
    % the per-class one: a class reaches counts up to its chain's total, and a
    % class declared empty and filled only by a switch (a ClassSwitch target, a
    % Cache hit/miss class) reaches them from a declared population of 0.
    % Bounding each class by its OWN population leaves reachable microstates
    % unenumerated -- every mixed one, when that population is 0 -- and the
    % check then passes vacuously on plainly order-dependent rates. csMatrix is
    % the class-switching mask link() has just built above; its connected
    % components are the chains. Without class switching each component is a
    % single class and Nvec falls back to pop.
    % a mask of a different width belongs to a different class set (a
    % fork-join transformation copies the model and then adds classes); ignore
    % it and fall back to per-class populations, as the chain code does.
    if ~isempty(self.csMatrix) && isequal(size(self.csMatrix), [K K])
        cs = self.csMatrix | self.csMatrix';
    else
        cs = logical(eye(K));
    end
    comp = zeros(1, K);
    ncomp = 0;
    for r = 1:K
        if comp(r) == 0
            ncomp = ncomp + 1;
            frontier = r;
            while ~isempty(frontier)
                v = frontier(1); frontier(1) = [];
                if comp(v) ~= 0, continue; end
                comp(v) = ncomp;
                frontier = [frontier, find(cs(v,:) & comp == 0)]; %#ok<AGROW>
            end
        end
    end
    Nvec = zeros(1, K);
    for r = 1:K
        Nvec(r) = sum(pop(comp == comp(r)));
    end
    % ntot caps the microstate LENGTH: with Nvec now chain-wide, summing it
    % inside checkPermInvariance would count each chain once per class and admit
    % microstates longer than the network can produce.
    ntot = sum(pop);
    for ind = 1:I
        nd = self.nodes{ind};
        if isa(nd, 'Queue') && SchedStrategy.toId(nd.schedStrategy) == SchedStrategy.OI ...
                && ~isempty(nd.svcRateFun)
            [ok, badc, partial] = nd.checkPermInvariance(Nvec, min(nd.cap, ntot));
            if ~ok
                line_error(mfilename, 'Order-independent (OI) station ''%s'' has a service rate function that is not permutation-invariant: mu(c) differs for a reordering of the microstate %s. Use SchedStrategy.PAS for order-dependent service, or disable this check with model.setChecks(false).', nd.getName(), mat2str(badc));
            end
            % (1): per-job rates must be non-negative. Runs second because
            % permutation invariance is what makes mu a function of the count
            % vector, which is what the increment test walks.
            [okm, badm, badr, partialm] = nd.checkRateMonotonicity(Nvec, min(nd.cap, ntot));
            partial = partial || partialm;
            if ~okm
                line_error(mfilename, 'Order-independent (OI) station ''%s'' has a service rate function with a negative per-job service rate: mu(c) DROPS when a class-%d job joins, at microstate %s. An OI rate must satisfy mu(c1..cj) >= mu(c1..c_{j-1}), so that every mu_j(c) is non-negative. A processor-sharing total rate (sum_j mu_{cj})/n has this shape whenever classes have different rates -- use SchedStrategy.PS for that station, or disable this check with model.setChecks(false).', nd.getName(), badr, mat2str(badm));
            elseif partial
                line_warning(mfilename, 'Order-independent (OI) station ''%s'': the permutation-invariance check was only partial because the reachable population is large; a subset of microstates was verified. To skip this check, call model.setChecks(false) before link().\n', nd.getName());
            end
        end
    end
end

%% DEBUGGING
if false
    % create java version of the network
    if isempty(self.obj)
        jnetwork = JLINE.from_line_network(self);
    else
        jnetwork = self.obj;
    end
    % compare sn data structures
    jsn = JLINE.from_jline_struct(jnetwork);
    fprintf(1,'* Comparison with Java NetworkStruct: ');
    bool = testJavaStruct(self.getName(),self.getStruct(),jsn);
end

end
