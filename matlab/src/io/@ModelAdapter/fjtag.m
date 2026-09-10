function [fjmodel, fjsn, fjclassmap, fjforkmap, fjjoinmap, fjbranchmap, fjtagmap] = fjtag(model)
% [FJMODEL, FJSN, FJCLASSMAP, FJFORKMAP, FJJOINMAP, FJBRANCHMAP, FJTAGMAP] = FJTAG(MODEL)
%
% Build a tag-augmented copy of a closed fork-join model for exact native
% analysis by SolverCTMC/SolverSSA. For each (fork f, class r) pair with
% matched join j, branch b=1..B and tag t=1..njobs(r), an auxiliary
% transient closed class A(f,r,b,t) with population 0 is created. The tag
% identifies the origin job: the fork firing (State.afterFJEvent) emits
% one sibling per branch in classes A(f,r,1..B,t) using the lowest free
% tag, and the join fires only when all siblings OF THE SAME TAG are
% buffered, releasing one class-r job (State.afterEventJoin). Identity
% matching is therefore exact even when siblings overtake each other
% across branches.
%
% Outputs:
% fjmodel      - augmented Network copy (forkStateful/isFJAugmented set)
% fjsn         - augmented NetworkStruct with sn.fjsync and nodeparam fj
%                blocks; run the CTMC/SSA algorithms on this struct
% fjclassmap(a) - original class of auxiliary class a (0 for originals)
% fjforkmap(a)  - fork node of auxiliary class a (0 for originals)
% fjjoinmap(a)  - join node of auxiliary class a (0 for originals)
% fjbranchmap(a)- branch index of auxiliary class a (0 for originals)
% fjtagmap(a)   - tag (origin-job slot) of auxiliary class a (0 for originals)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

sn = model.getStruct;
sn_fj_validate(sn);

K = sn.nclasses;
I = sn.nnodes;
Vnodes = cellsum(sn.nodevisits);

fjmodel = model.copy();
fjmodel.allowReplace = true;
P = fjmodel.getLinkedRoutingMatrix;
if isempty(P)
    line_error(mfilename,'The native CTMC/SSA fork-join implementation requires the routing topology to be generated using Network.link.');
end
fjmodel.resetNetwork(true);
fjmodel.resetStruct();
fjmodel.forkStateful = true;
fjmodel.isFJAugmented = true;
% the linked routing matrix is indexed over user nodes only (resetNetwork
% deletes auto-inserted ClassSwitch nodes, relink re-creates them), while
% sn indices include them: keep the two spaces distinct
Ip = size(P{1,1},1);

fjclassmap = zeros(1,K);
fjforkmap = zeros(1,K);
fjjoinmap = zeros(1,K);
fjbranchmap = zeros(1,K);
fjtagmap = zeros(1,K);

forkIndexes = find(sn.nodetype == NodeType.Fork)';
% per-fork bookkeeping for the sn post-edits
forkinfo = {}; % rows: {f, j, r, branchheads, branchsets, auxmatrix (BxT)}

for f=forkIndexes
    j = find(sn.fj(f,:));
    % tasks emitted per output link at each fork firing. `w` is the node-wide
    % count and `wlink` the per-destination one; they differ only when the fork
    % declares a variable forking level, and sn_fj_validate has already refused
    % the cases this construction cannot carry (a random degree, a non-integer
    % count, an uncertain branch under a standard Join).
    w = 1;
    if isfield(sn.nodeparam{f},'fanOut') && ~isempty(sn.nodeparam{f}.fanOut)
        w = round(sn.nodeparam{f}.fanOut(1));
    end
    wlink = [];
    if isfield(sn.nodeparam{f},'fanOutLink') && ~isempty(sn.nodeparam{f}.fanOutLink)
        wlink = round(sn.nodeparam{f}.fanOutLink);
    end
    for r=find(Vnodes(f,:)>0)
        % branch heads: nodes receiving class r directly from the fork
        branchheads = [];
        for jnd=1:I
            if sn.rtnodes((f-1)*K+r,(jnd-1)*K+r) > 0
                branchheads(end+1) = jnd; %#ok<AGROW>
            end
        end
        B = length(branchheads);
        if B < 2
            line_error(mfilename,'Degenerate forks with a single output link are not supported by the native CTMC/SSA fork-join implementation.');
        end
        % branch discovery: class-r BFS closure from each head up to the join
        branchsets = cell(1,B);
        for b=1:B
            visitset = branchheads(b);
            frontier = branchheads(b);
            while ~isempty(frontier)
                cn = frontier(1); frontier(1) = [];
                if sn.nodetype(cn) == NodeType.Fork
                    line_error(mfilename,'Nested fork-join is not supported by the native CTMC/SSA fork-join implementation.');
                end
                if sn.nodetype(cn) == NodeType.Join && cn ~= j
                    line_error(mfilename,'Overlapping fork-join pairs are not supported by the native CTMC/SSA fork-join implementation.');
                end
                % see _kb/12-interfaces-and-docs.md (@ModelAdapter: JMT export helpers) for rationale
                if cn > Ip
                    line_error(mfilename,'Class switching between fork and join is not supported by the native CTMC/SSA fork-join implementation.');
                end
                for jnd=1:I
                    for s=1:K
                        if sn.rtnodes((cn-1)*K+r,(jnd-1)*K+s) > 0
                            if s ~= r
                                line_error(mfilename,'Class switching between fork and join is not supported by the native CTMC/SSA fork-join implementation.');
                            end
                            if jnd ~= j && ~ismember(jnd, visitset)
                                visitset(end+1) = jnd; %#ok<AGROW>
                                frontier(end+1) = jnd; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
            % trap check: every branch node must reach the join within the branch
            canreach = j;
            changed = true;
            while changed
                changed = false;
                for cn=visitset
                    if ~ismember(cn, canreach)
                        for jnd=canreach
                            if sn.rtnodes((cn-1)*K+r,(jnd-1)*K+r) > 0
                                canreach(end+1) = cn; %#ok<AGROW>
                                changed = true;
                                break
                            end
                        end
                    end
                end
            end
            if ~all(ismember(visitset, canreach))
                line_error(mfilename,'Fork branches from which the Join is unreachable are not supported by the native CTMC/SSA fork-join implementation.');
            end
            branchsets{b} = visitset;
        end
        % tag pool size = maximum number of concurrently outstanding forked
        % jobs of class r = population of its chain (class switching outside
        % the fork-join section can concentrate the whole chain population
        % in class r, e.g. a class switch on the edge into the fork)
        c = find(sn.chains(:,r), 1);
        T = round(sum(sn.njobs(sn.chains(c,:))));
        if ~isfinite(T)
            line_error(mfilename,'Chains with infinite population routed through a Fork are not supported by the native CTMC/SSA fork-join implementation.');
        end
        auxmatrix = zeros(B,T);
        for t=1:T
            for b=1:B
                auxname = sprintf('%s_f%d_b%d_t%d', sn.classnames{r}, f, b, t);
                auxclass = ClosedClass(fjmodel, auxname, 0, fjmodel.stations{sn.refstat(r)}, sn.classprio(r));
                a = auxclass.index;
                auxmatrix(b,t) = a;
                fjclassmap(a) = r;
                fjforkmap(a) = f;
                fjjoinmap(a) = j;
                fjbranchmap(a) = b;
                fjtagmap(a) = t;
                % sibling service on the branch copies the original class
                for cn=branchsets{b}
                    if sn.isstation(cn) && sn.nodetype(cn) ~= NodeType.Join
                        svc = model.nodes{cn}.getService(model.classes{r});
                        fjmodel.nodes{cn}.setService(auxclass, svc.copy());
                    end
                end
                % register the auxiliary class at the join input section
                fjmodel.nodes{j}.setStrategy(auxclass, JoinStrategy.STD);
                fjmodel.nodes{j}.setRequired(auxclass, -1);
                % sibling routing: copy the class-r branch routing; the
                % auxiliary class terminates at the join (no outgoing row)
                P{auxclass,auxclass} = zeros(Ip);
                for cn=branchsets{b}
                    P{auxclass,auxclass}(cn,:) = P{r,r}(cn,:);
                end
            end
        end
        forkinfo(end+1,:) = {f, j, r, branchheads, branchsets, auxmatrix, w, wlink}; %#ok<AGROW>
    end
end

% complete the routing cell array over all class pairs
Kaug = length(fjmodel.classes);
for a=1:Kaug
    for b=1:Kaug
        if size(P,1)<a || size(P,2)<b || isempty(P{a,b})
            P{a,b} = zeros(Ip);
        end
    end
end
% re-initialize output dispatchers of non-station nodes over the enlarged
% class set (resetNetwork only re-initializes station dispatchers); the
% fork out-links are restored by relink from the unchanged P{r,r} rows
for nd = 1:length(fjmodel.nodes)
    if ~isa(fjmodel.nodes{nd}, 'Station') && isprop(fjmodel.nodes{nd}, 'output') ...
            && ~isempty(fjmodel.nodes{nd}.output) && ismethod(fjmodel.nodes{nd}.output, 'initDispatcherJobClasses')
        fjmodel.nodes{nd}.output.initDispatcherJobClasses(fjmodel.classes);
    end
end
fjmodel.relink(P);

% see _kb/12-interfaces-and-docs.md (@ModelAdapter: JMT export helpers) for rationale
csm = fjmodel.csMatrix;
if isempty(csm)
    csm = logical(eye(Kaug));
end
if size(csm,1) < Kaug
    csm(Kaug,Kaug) = false;
end
for a=1:Kaug
    csm(a,a) = true;
    if fjclassmap(a) > 0
        csm(fjclassmap(a),a) = true;
        csm(a,fjclassmap(a)) = true;
    end
end
fjmodel.csMatrix = logical(csm);

% re-initialize the default state: the copy inherits the original model's
% initialization flag and states, which have pre-augmentation class widths
fjmodel.initDefault();

fjsn = fjmodel.getStruct();

% ---- sn post-edits ----

% see _kb/12-interfaces-and-docs.md (@ModelAdapter: JMT export helpers) for rationale
for r=1:K
    corig = find(sn.chains(:,r), 1);
    cnew = find(fjsn.chains(:,r), 1);
    for isfnew=1:fjsn.nstateful
        ind = fjsn.statefulToNode(isfnew);
        if sn.isstateful(ind)
            fjsn.visits{cnew}(isfnew,r) = sn.visits{corig}(sn.nodeToStateful(ind),r);
        else % stateful Fork: nonzero marker for capacity gating
            fjsn.visits{cnew}(isfnew,r) = double(Vnodes(ind,r) > 0);
        end
    end
    nshared = min(size(fjsn.nodevisits{cnew},1), size(sn.nodevisits{corig},1));
    fjsn.nodevisits{cnew}(1:nshared,r) = sn.nodevisits{corig}(1:nshared,r);
end

% auxiliary-class visits: engines use them only as a zero-versus-nonzero
% capacity gate, so set 1 on the branch support set and 0 elsewhere
for row=1:size(forkinfo,1)
    [f, j, r, branchheads, branchsets, auxmatrix, w, wlink] = forkinfo{row,:};
    B = size(auxmatrix,1);
    T = size(auxmatrix,2);
    cnew = find(fjsn.chains(:,r), 1);
    for b=1:B
        support = [branchsets{b}, j];
        for t=1:T
            a = auxmatrix(b,t);
            fjsn.visits{cnew}(:,a) = 0;
            fjsn.nodevisits{cnew}(:,a) = 0;
            for cn=support
                fjsn.nodevisits{cnew}(cn,a) = 1;
                if fjsn.isstateful(cn)
                    fjsn.visits{cnew}(fjsn.nodeToStateful(cn),a) = 1;
                end
            end
            % Capacity: each auxiliary class holds at most the tasks THIS
            % BRANCH is sent, network-wide (STD join, one tag at a time). It is
            % the branch's own count and not the node-wide one, because a fork
            % that sends 1 down one link and 3 down another would otherwise cap
            % the second branch at the mean and deadlock the chain.
            if isempty(wlink)
                wb = w;
            else
                wb = wlink(branchheads(b), r);
            end
            fjsn.classcap(:,a) = 0;
            for cn=support
                if fjsn.isstation(cn)
                    fjsn.classcap(fjsn.nodeToStation(cn),a) = wb;
                end
            end
        end
    end
    % nodeparam fj blocks read by State.afterEventJoin/afterFJEvent
    if ~isstruct(fjsn.nodeparam{f})
        fjsn.nodeparam{f} = struct();
    end
    if ~isfield(fjsn.nodeparam{f},'fj')
        fjsn.nodeparam{f}.fj = struct('classes',[],'joins',[],'auxmatrix',{cell(1,Kaug)},'branchheads',{cell(1,Kaug)});
    end
    fjsn.nodeparam{f}.fj.classes(end+1) = r; % original classes forked here
    fjsn.nodeparam{f}.fj.joins(end+1) = j;
    fjsn.nodeparam{f}.fj.auxmatrix{r} = auxmatrix;
    fjsn.nodeparam{f}.fj.branchheads{r} = branchheads;
    if ~isstruct(fjsn.nodeparam{j})
        fjsn.nodeparam{j} = struct();
    end
    if ~isfield(fjsn.nodeparam{j},'fj')
        fjsn.nodeparam{j}.fj = struct('fork',f,'origclasses',[],'auxmatrix',{cell(1,Kaug)},'required',{cell(1,Kaug)});
    end
    fjsn.nodeparam{j}.fj.origclasses(end+1) = r;
    fjsn.nodeparam{j}.fj.auxmatrix{r} = auxmatrix;
    % Siblings of branch b a firing consumes. Under STD it is the tasks that
    % branch was sent, which is the per-destination count when the fork declares
    % one and the node-wide count otherwise. Under PARTIAL the Join's own fanIn
    % is what it waits for, so the quorum lowers this count rather than
    % replacing the mechanism.
    reqb = w*ones(B,1);
    if ~isempty(wlink)
        for bb=1:B
            reqb(bb) = wlink(branchheads(bb), r);
        end
    end
    if isfield(sn.nodeparam{j},'joinStrategy') && length(sn.nodeparam{j}.joinStrategy) >= r && ...
            ~isempty(sn.nodeparam{j}.joinStrategy{r}) && sn.nodeparam{j}.joinStrategy{r} == JoinStrategy.PARTIAL && ...
            isfield(sn.nodeparam{j},'fanIn') && length(sn.nodeparam{j}.fanIn) >= r && ~isempty(sn.nodeparam{j}.fanIn{r})
        fq = sn.nodeparam{j}.fanIn{r};
        if isscalar(fq)
            fq = fq*ones(B,1);
        end
        % both as COLUMNS: a Bx1 against a 1xB broadcasts to BxB in MATLAB,
        % which is not a per-branch requirement but a matrix nothing can read
        fq = fq(:);
        reqb = min(reqb(:), fq(1:B));
    end
    fjsn.nodeparam{j}.fj.required{r} = reqb;
end

% fork firing synchronizations: one entry per (fork, class, tag)
fjsn.fjsync = {};
for row=1:size(forkinfo,1)
    [f, j, r, branchheads, ~, auxmatrix, w, wlink] = forkinfo{row,:};
    T = size(auxmatrix,2);
    for t=1:T
        entry = struct();
        entry.active{1} = Event(EventType.FIRE, f, r);
        entry.fork = f;
        entry.join = j;
        entry.class = r;
        entry.tag = t;
        entry.branchheads = branchheads;
        entry.auxclasses = auxmatrix(:,t)';
        entry.auxall = auxmatrix; % B x T, for the tag-occupancy scan
        % Siblings emitted per branch, per destination when the fork sends
        % different counts down different links. It stays the SCALAR whenever
        % every branch agrees, so a plain fork takes exactly the code path it
        % always took.
        if isempty(wlink)
            entry.weight = w;
        else
            wvec = arrayfun(@(bb) wlink(bb, r), branchheads(:)');
            if all(wvec == wvec(1))
                entry.weight = wvec(1);
            else
                entry.weight = wvec;
            end
        end
        % branch activation probability, one per branch; all ones on a fork whose
        % branches are certain, which is every fork the CTMC/SSA path accepts
        % under a standard Join
        if isfield(sn.nodeparam{f},'fanOutProb') && ~isempty(sn.nodeparam{f}.fanOutProb)
            pvec = arrayfun(@(bb) sn.nodeparam{f}.fanOutProb(bb, r), branchheads(:)');
            if all(pvec == 1)
                entry.prob = 1.0;
            else
                entry.prob = pvec;
            end
        else
            entry.prob = 1.0;
        end
        fjsn.fjsync{end+1,1} = entry;
    end
end

% auxiliary class map, consumed by capacity gating (e.g. the SSA preamble
% must not cap sibling multiplicities by the chain job population)
fjsn.fjclassmap = fjclassmap;

fjmodel.sn = fjsn;

end
