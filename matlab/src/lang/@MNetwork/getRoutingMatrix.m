function [rt,rtnodes,conn,chains,rtNodesByClass,rtNodesByStation] = getRoutingMatrix(self, arvRates)
% [RT,RTNODES,CONNMATRIX,CHAINS,RTNODESBYCLASS,RTNODESBYSTATION] = GETROUTINGMATRIX(ARVRATES)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if self.hasStruct
    sn = self.sn;
    idxSource = find(sn.nodetype == NodeType.Source);
    idxSink = find(sn.nodetype == NodeType.Sink);
    indexOpenClasses = find(isinf(sn.njobs));
    hasOpen = ~isempty(indexOpenClasses);
    if nargin<2
        arvRates = sn.rates(sn.nodeToStation(idxSource),:);
    end
    nodes = self.nodes;
    conn = sn.connmatrix;
    I = sn.nnodes;
    K = sn.nclasses;
    NK = sn.njobs;
    statefulNodes = find(sn.isstateful)';
else
    sn = self.sn;
    idxSource = self.getIndexSourceNode;
    idxSink = self.getIndexSinkNode;
    indexOpenClasses = self.getIndexOpenClasses;
    nodes = self.nodes;

    conn = self.getConnectionMatrix;
    sn.connmatrix = conn;
    hasOpen = hasOpenClasses(self);

    I = self.getNumberOfNodes;
    K = self.getNumberOfClasses;
    NK = self.getNumberOfJobs;
    for ind=1:I
        for k=1:K
            if isempty(self.nodes{ind}.output.outputStrategy{k})
                sn.routing(ind,k) = RoutingStrategy.DISABLED;
            else
                sn.routing(ind,k) = RoutingStrategy.fromText(nodes{ind}.output.outputStrategy{k}{2});
            end
        end
    end
    statefulNodes = self.getIndexStatefulNodes;
end

rtnodes = zeros(I*K);
rtNodesByClass = {};
rtNodesByStation = {};
% The first loop considers the class at which a job enters the
% target station
for ind=1:I
    node_i = nodes{ind};
    outputStrategy = node_i.output.outputStrategy;
    switch class(node_i.output)
        case 'Forker'
            % Fork nodes send to ALL branches with probability 1.0 each
            % (not 1.0/fanout). A Fork with arrival rate lambda sends
            % lambda to EACH outgoing branch, not lambda/fanout.
            for k=1:K
                outputStrategy_k = outputStrategy{k};
                switch sn.routing(ind,k)
                    case RoutingStrategy.PROB
                        % Use output strategy destinations (includes ClassSwitch nodes)
                        if ~isempty(outputStrategy_k) && length(outputStrategy_k) >= 3 && ~isempty(outputStrategy_k{end})
                            for t=1:length(outputStrategy_k{end})
                                jnd = outputStrategy_k{end}{t}{1}.index;
                                rtnodes((ind-1)*K+k,(jnd-1)*K+k) = 1.0;
                            end
                        else
                            % Fallback: uniform to all connected nodes
                            for jnd=1:I
                                if conn(ind,jnd)>0
                                    rtnodes((ind-1)*K+k,(jnd-1)*K+k) = 1.0;
                                end
                            end
                        end
                    otherwise
                        % Non-PROB routing: route to all connected nodes
                        for jnd=1:I
                            if conn(ind,jnd)>0
                                rtnodes((ind-1)*K+k,(jnd-1)*K+k) = 1.0;
                            end
                        end
                end
            end
        otherwise
            isSink_i = isa(node_i,'Sink');
            for k=1:K
                outputStrategy_k = outputStrategy{k};
                switch sn.routing(ind,k)
                    case RoutingStrategy.PROB
                        if isinf(NK(k)) || ~isSink_i
                            for t=1:length(outputStrategy_k{end}) % for all outgoing links
                                jnd = outputStrategy_k{end}{t}{1}.index;
                                rtnodes((ind-1)*K+k,(jnd-1)*K+k) = outputStrategy_k{end}{t}{2};
                            end
                        end
                    case RoutingStrategy.DISABLED
                        for jnd=1:I
                            if conn(ind,jnd)>0
                                if isa(self.classes{k},'SelfLoopingClass')
                                    rtnodes((ind-1)*K+k,(jnd-1)*K+k) = 0;
                                elseif isa(self.classes{k},'Signal') || isa(self.classes{k},'OpenSignal') || isa(self.classes{k},'ClosedSignal')
                                    % Signal classes with DISABLED routing should not route
                                    % (their routing is set explicitly via the P matrix)
                                    rtnodes((ind-1)*K+k,(jnd-1)*K+k) = 0;
                                else
                                    % we set this to be non-zero as otherwise the
                                    % classes that do not visit a classswitch are
                                    % misconfigured in JMT
                                    rtnodes((ind-1)*K+k,(jnd-1)*K+k) = 1/sum(conn(ind,:));
                                end
                            end
                        end
                    case {RoutingStrategy.RAND, RoutingStrategy.RROBIN, RoutingStrategy.JSQ}
                        if isinf(NK(k)) % open class
                            for jnd=1:I
                                if conn(ind,jnd)>0
                                    rtnodes((ind-1)*K+k,(jnd-1)*K+k)=1/sum(conn(ind,:));
                                end
                            end
                        elseif ~isa(node_i,'Source') && ~isSink_i % don't route closed classes out of source nodes
                            connectionsClosed = conn;
                            if connectionsClosed(ind,idxSink)
                                connectionsClosed(ind,idxSink) = 0;
                            end
                            for jnd=1:I
                                if connectionsClosed(ind,jnd)>0
                                    rtnodes((ind-1)*K+k,(jnd-1)*K+k)=1/(sum(connectionsClosed(ind,:)));
                                end
                            end
                        end
                    case RoutingStrategy.WRROBIN
                        % Weighted round-robin: compute routing probabilities from weights
                        % Weights are stored in outputStrategy_k{3} as {destination, weight} pairs
                        if ~isempty(outputStrategy_k) && length(outputStrategy_k) >= 3 && ~isempty(outputStrategy_k{3})
                            % Compute total weight (excluding sink for closed classes)
                            totalWeight = 0;
                            for t=1:length(outputStrategy_k{3})
                                dest = outputStrategy_k{3}{t}{1};
                                destIdx = dest.index;
                                weight = outputStrategy_k{3}{t}{2};
                                % For closed classes, exclude sink
                                if ~isinf(NK(k)) && ~isempty(idxSink) && destIdx == idxSink
                                    continue;
                                end
                                totalWeight = totalWeight + weight;
                            end
                            % Set routing probabilities based on weights
                            if totalWeight > 0
                                for t=1:length(outputStrategy_k{3})
                                    dest = outputStrategy_k{3}{t}{1};
                                    destIdx = dest.index;
                                    weight = outputStrategy_k{3}{t}{2};
                                    % For closed classes, exclude sink
                                    if ~isinf(NK(k)) && ~isempty(idxSink) && destIdx == idxSink
                                        continue;
                                    end
                                    rtnodes((ind-1)*K+k,(destIdx-1)*K+k) = weight / totalWeight;
                                end
                            else
                                % Fallback to uniform if no weights sum to positive
                                if isinf(NK(k))
                                    for jnd=1:I
                                        if conn(ind,jnd)>0
                                            rtnodes((ind-1)*K+k,(jnd-1)*K+k)=1/sum(conn(ind,:));
                                        end
                                    end
                                elseif ~isa(node_i,'Source') && ~isSink_i
                                    connectionsClosed = conn;
                                    if ~isempty(idxSink) && connectionsClosed(ind,idxSink)
                                        connectionsClosed(ind,idxSink) = 0;
                                    end
                                    for jnd=1:I
                                        if connectionsClosed(ind,jnd)>0
                                            rtnodes((ind-1)*K+k,(jnd-1)*K+k)=1/(sum(connectionsClosed(ind,:)));
                                        end
                                    end
                                end
                            end
                        else
                            % Fallback to uniform if no weight data available
                            if isinf(NK(k))
                                for jnd=1:I
                                    if conn(ind,jnd)>0
                                        rtnodes((ind-1)*K+k,(jnd-1)*K+k)=1/sum(conn(ind,:));
                                    end
                                end
                            elseif ~isa(node_i,'Source') && ~isSink_i
                                connectionsClosed = conn;
                                if ~isempty(idxSink) && connectionsClosed(ind,idxSink)
                                    connectionsClosed(ind,idxSink) = 0;
                                end
                                for jnd=1:I
                                    if connectionsClosed(ind,jnd)>0
                                        rtnodes((ind-1)*K+k,(jnd-1)*K+k)=1/(sum(connectionsClosed(ind,:)));
                                    end
                                end
                            end
                        end
                    otherwise
                        for jnd=1:I
                            if conn(ind,jnd)>0
                                rtnodes((ind-1)*K+k,(jnd-1)*K+k) = GlobalConstants.FineTol;
                            end
                        end
                        %line_error(mfilename,[outputStrategy_k{2},' routing policy is not yet supported.']);
                end
            end
    end
end

% The second loop corrects the first one at nodes that change
% the class of the job in the service section.
for ind=1:I % source
    if isa(nodes{ind}.server,'StatelessClassSwitcher')
        Pi = rtnodes(((ind-1)*K+1):ind*K,:);
        Pcs = nodes{ind}.server.csFun(1:K,1:K);
        rtnodes(((ind-1)*K+1):ind*K,:) = 0;
        for jnd=1:I % destination
            Pij = Pi(:,((jnd-1)*K+1):jnd*K); %Pij(r,s)
            % Find the routing probability section determined by the router section in the first loop
            rtnodes(((ind-1)*K+1) : ((ind-1)*K+K),(jnd-1)*K+(1:K)) = Pcs.*repmat(diag(Pij)',K,1);
        end
    elseif isa(nodes{ind}.server,'StatefulClassSwitcher')
        Pi = rtnodes(((ind-1)*K+1):ind*K,:);
        for r=1:K
            for s=1:K
                Pcs(r,s) = nodes{ind}.server.csFun(r,s,[],[]); % get csmask
            end
        end
        rtnodes(((ind-1)*K+1):ind*K,:) = 0;
        if isa(nodes{ind}.server,'CacheClassSwitcher')
            for r=1:K
                if (isempty(find(r == nodes{ind}.server.hitClass)) && isempty(find(r == nodes{ind}.server.missClass)))
                    Pcs(r,:) = Pcs(r,:)/sum(Pcs(r,:));
                end
            end

            for r=1:K
                if (isempty(find(r == nodes{ind}.server.hitClass)) && isempty(find(r == nodes{ind}.server.missClass)))
                    for jnd=1:I % destination
                        for s=1:K
                            Pi((ind-1)*K+r,(jnd-1)*K+s) = 0;
                        end
                    end
                end
            end

            for r=1:K
                if length(nodes{ind}.server.actualHitProb)>=r && length(nodes{ind}.server.hitClass)>=r
                    ph = nodes{ind}.server.actualHitProb(r);
                    pm = nodes{ind}.server.actualMissProb(r);
                    h = nodes{ind}.server.hitClass(r);
                    m = nodes{ind}.server.missClass(r);
                    rtnodes((ind-1)*K+r,(ind-1)*K+h) = ph;
                    rtnodes((ind-1)*K+r,(ind-1)*K+m) = pm;
                else
                    if length(nodes{ind}.server.hitClass)>=r
                        h = nodes{ind}.server.hitClass(r);
                        m = nodes{ind}.server.missClass(r);
                        if h > 0 && m > 0  % Only set NaN if hitClass/missClass are defined for this class
                            rtnodes((ind-1)*K+r,(ind-1)*K+h) = NaN;
                            rtnodes((ind-1)*K+r,(ind-1)*K+m) = NaN;
                        end
                    end
                end
            end
            for jnd=1:I % destination
                Pij = Pi(1:K,((jnd-1)*K+1):jnd*K); %Pij(r,s)
                for r=1:K
                    if ~(isempty(find(r == nodes{ind}.server.hitClass)) && isempty(find(r == nodes{ind}.server.missClass)))
                        for s=1:K
                            % Find the routing probability section determined by the router section in the first loop
                            %Pnodes(((i-1)*K+1):i*K,((j-1)*K+1):j*K) = Pcs*Pij;
                            rtnodes((ind-1)*K+r,(jnd-1)*K+s) = Pcs(r,s)*Pij(s,s);
                        end
                    end
                end
            end
        end
    end
end

% Handle self-looping classes
for r=1:K
    if isa(self.classes{r},'SelfLoopingClass')
        sl_refstat = self.classes{r}.refstat.index;
        for ind=1:I % source
            for jnd=1:I % destination
                for s=1:K
                    rtnodes((ind-1)*K+r,(jnd-1)*K+s) = 0; % zero any prior setting
                end
            end
        end
        for s=1:K
            if s==r
                rtnodes((sl_refstat-1)*K+r,(sl_refstat-1)*K+s) = 1;
            end
        end
    end
end

% ignore all chains containing a Pnodes column that sums to 0,
% since these are classes that cannot arrive to the node
% unless this column belongs to the source
colsToIgnore = find(sum(rtnodes,1)==0);
if hasOpen
    colsToIgnore = setdiff(colsToIgnore,(idxSource-1)*K+(1:K));
end

% We route back from the sink to the source. Since open classes
% have an infinite population, if there is a class switch QN
% with the following chains
% Source -> (A or B) -> C -> Sink
% Source -> D -> Sink
% We can re-route class C into the source either as A or B or C.
% We here re-route back as C and leave for the chain analyzer
% to detect that C is in a chain with A and B and change this
% part.

csmask = self.csMatrix;
if isempty(csmask) % models not built with link(P)
    % rtnodes+rtnodes' makes the matrix undirected, it helps grouping
    % transient states with all the recurrent states their can reach
    % into a single "chain" of the LINE model
    [C,WCC]=weaklyconncomp(rtnodes+rtnodes');
    WCC(colsToIgnore) = 0;
    chainCandidates = cell(1,C);
    for r=1:C
        chainCandidates{r} = find(WCC==r);
    end

    chains = false(length(chainCandidates),K); % columns are classes, rows are chains
    tmax = 0;
    for t=1:length(chainCandidates)
        if length(chainCandidates{t})>1
            tmax = tmax+1;
            chains(tmax,(mod(chainCandidates{t}-1,K)+1))=true;
        end
    end
    chains(tmax+1:end,:)=[];
else % this is faster since csmask is smaller than rtnodes
    [C,WCC] = weaklyconncomp(csmask+csmask');

    chainCandidates = cell(1,C);
    for c=1:C
        chainCandidates{c} = find(WCC==c);
    end

    chains = false(length(chainCandidates),K);
    for t=1:length(chainCandidates)
        chains(t,chainCandidates{t}) = true;
    end
end

chains = unique(chains,'rows');
try
    chains = sortrows(chains,'descend');
catch % old MATLABs
    chains = sortrows(chains);
end

% overlapping chains can occur if a class that enters in a Cache
% reappars downhill in the network in rtnodes in that case we merge
% the two chains
splitChains = [];
for col=find(sum(chains,1)>1)
    rows = find(chains(:,col));
    for r=rows(2:end)'
        chains(rows(1),:) = chains(rows(1),:) | chains(r,:);
        splitChains(end+1)=r;
    end
end
chains(splitChains,:)=[];

% We now obtain the routing matrix P by ignoring the non-stateful
% nodes and calculating by the stochastic complement method the
% correct transition probabilities, that includes the effects
% of the non-stateful nodes (e.g., ClassSwitch)
statefulNodesClasses = [];
for ind=statefulNodes %#ok<FXSET>
    statefulNodesClasses(end+1:end+K)= ((ind-1)*K+1):(ind*K);
end

% this routes open classes back from the sink into the source
% it will not work with non-renewal arrivals as it choses in which open
% class to reroute a job with probability depending on the arrival rates
if hasOpen
    arvRates(isnan(arvRates)) = 0;
    for s=indexOpenClasses
        s_chain = find(chains(:,s));
        others_in_chain = find(chains(s_chain,:)); %#ok<FNDSB>
        rtnodes((idxSink-1)*K+others_in_chain,(idxSource-1)*K+others_in_chain) = repmat(arvRates(others_in_chain)/sum(arvRates(others_in_chain)),length(others_in_chain),1);
    end
end

%% Hide the nodes that are not stateful
rt = dtmc_stochcomp(rtnodes,statefulNodesClasses);

% verify irreducibility
% disabled as it casts warning in many models that are seemingly fine
%try
%    eigen = eig(rt); % exception if rt has NaN or Inf
%    if sum(eigen>=1)>1
%        line_warning(mfilename, 'Solutions may be invalid, the routing matrix is reducible. The path of two or more classes form independent non-interacting models.');
%    end
%end
sn.rt = rt;
sn.chains = chains;
self.sn = sn;

%% Compute optional outputs
if nargout >= 5
    rtNodesByClass = cellzeros(K,K,I,I);
    for ind=1:I
        for jnd=1:I
            for r=1:K
                for s=1:K
                    rtNodesByClass{s,r}(ind,jnd) = rtnodes((ind-1)*K+s,(jnd-1)*K+r);
                end
            end
        end
    end
end

if nargout >= 6
    rtNodesByStation = cellzeros(I,I,K,K);
    for ind=1:I %#ok<FXSET>
        for jnd=1:I
            for r=1:K
                for s=1:K
                    rtNodesByStation{ind,jnd}(r,s) = rtnodes((ind-1)*K+s,(jnd-1)*K+r);
                end
            end
        end
    end
end
end