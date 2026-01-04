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
    % in this case it is possible that P is linear but just because the
    % routing is state-dependent and therefore some zero entries are
    % actually unspecified
    %cacheNodes = find(cellfun(@(c) isa(c,'Cache'), self.getStatefulNodes));
    for ind=1:I
        switch class(self.nodes{ind})
            case 'Cache'
                % note that since a cache needs to distinguish hits and
                % misses, it needs to do class-switch unless the model is
                % degenerate
                isLinearP = false;
                if self.nodes{ind}.server.hitClass == self.nodes{ind}.server.missClass
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
% eye(K) is because any job that travels to an autoAdded class
% switch stays in the same class prior to reaching the ClassSwitch
% and anyway the diagonal of csMatrix is irrelevant for the chains
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
                    end
                    P{s,s}(csid(ind,jnd),jnd) = 1;
                end
            end
        end
    end
end

connected = zeros(Ip);
nodes = self.nodes;
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

% check if the probability out of any node sums to >1.0
pSum = cellsum(P);
isAboveOne = pSum > 1.0 + GlobalConstants.FineTol;
if any(isAboveOne)
    for ind=find(isAboveOne)
        if SchedStrategy.toId(self.nodes{ind}.schedStrategy) ~= SchedStrategy.FORK
            line_error(mfilename,sprintf('The total routing probability for jobs leaving node %s in class %s is greater than 1.0.',self.nodes{ind}.name,self.classes{r}.name));
        end
        %        elseif pSum < 1.0 - GlobalConstants.FineTol % we cannot check this case as class r may not reach station i, in which case its outgoing routing prob is zero
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
