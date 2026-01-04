function [Q,stateSpace,stateSpaceAggr,Dfilt,arvRates,depRates,sn]=solver_ctmc(sn,options)
% [Q,SS,SSQ,DFILT,ARVRATES,DEPRATES,QN]=SOLVER_CTMC(QN,OPTIONS)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% generate state space
%nnodes = sn.nnodes;
nstateful = sn.nstateful;
nclasses = sn.nclasses;
sync = sn.sync;
A = length(sync);
csmask = sn.csmask;

line_debug('CTMC solver starting: nstateful=%d, nclasses=%d, sync_events=%d', nstateful, nclasses, A);

if ~isfield(options.config, 'hide_immediate')
    options.config.hide_immediate = true;
end

if ~isfield(options.config, 'state_space_gen')
    options.config.state_space_gen = 'default';
end

%% generate state spaces, detailed and aggregate
switch options.config.state_space_gen
    case 'reachable' % does not handle open models yet (no cutoff)
        line_debug('Using reachable state space generation, calling ctmc_ssg_reachability');
        [stateSpace, stateSpaceAggr, stateSpaceHashed,~,sn] = ctmc_ssg_reachability(sn,options);
    case {'default','full'}
        line_debug('Using full state space generation, calling ctmc_ssg');
        [stateSpace, stateSpaceAggr, stateSpaceHashed,~,sn] = ctmc_ssg(sn,options);
end

line_debug('State space generated: %d states', size(stateSpaceHashed,1));

%%
Q = sparse(eye(size(stateSpaceHashed,1))); % the diagonal elements will be removed later
Dfilt = cell(1,A);
for a=1:A
    Dfilt{a} = 0*Q;
end
local = sn.nnodes+1; % passive action

% SPN code
% Adj_t = zeros(size(SSh,1),size(SSh,1));
% Adj_m = zeros(size(SSh,1),size(SSh,1));
% if ~isempty(Adj) && ~isempty(ST)
%    edges = adj_to_mat(Adj);
% end

%% for all synchronizations
for a=1:A
    stateCell = cell(nstateful,1);
    %sn.sync{a}.active{1}.print
    for s=1:size(stateSpaceHashed,1)
        %[a,s]
        state = stateSpaceHashed(s,:);
        % SPN code
        %         ustate = stateSpace(s,:);
        %         state_pn = [];
        %         for st=1:length(ustate)
        %             if ~isempty(sn.varsparam{st}) && isfield(sn.varsparam{st}, 'nodeToPlace')
        %                 state_pn(sn.varsparam{st}.nodeToPlace) = ustate(st);
        %             end
        %         end

        % update state cell array and SSq
        for ind = 1:sn.nnodes
            if sn.isstateful(ind)
                isf = sn.nodeToStateful(ind);
                stateCell{isf} = sn.space{isf}(state(isf),:);
                %                if sn.isstation(ind)
                %                    ist = sn.nodeToStation(ind);
                %                    [~,nir] = State.toMarginal(sn,ind,stateCell{isf});
                %                end
            end
        end
        node_a = sync{a}.active{1}.node;
        state_a = state(sn.nodeToStateful(node_a));
        class_a = sync{a}.active{1}.class;
        event_a = sync{a}.active{1}.event;
        [new_state_a, rate_a] = State.afterEventHashed( sn, node_a, state_a, event_a, class_a);
        % SPN code:
        %[new_state_a, rate_a,~,trans_a, modes_a] = State.afterEventHashed( qn, node_a, state_a, event_a, class_a);

        %% debugging block
        %        if true%options.verbose == 2
        %            line_printf('---\n');
        %            sync{a}.active{1}.print,
        %        end
        %%
        if new_state_a == -1 % hash not found
            continue
        end
        for ia=1:length(new_state_a)
            if rate_a(ia)>0
                % SPN code:
                %if rate_a(ia)>0 || modes_a(ia) > 0
                node_p = sync{a}.passive{1}.node;
                if node_p ~= local
                    % Skip if the active transition hash was not found
                    if new_state_a(ia) == -1
                        continue
                    end
                    state_p = state(sn.nodeToStateful(node_p));
                    class_p = sync{a}.passive{1}.class;
                    event_p = sync{a}.passive{1}.event;

                    % SPN code:
                    %                     enabled = 0;
                    %                     if ia <= length(trans_a)
                    %                         % check if other input places of the transition contains as many token as the multiplicity of the input arcs
                    %                         tr = trans_a(ia);
                    %                         mode = modes_a(ia);
                    %                         bmatrix = sn.varsparam{tr}.back(:,mode);
                    %                         inmatrix = sn.varsparam{tr}.inh(:,mode);
                    %                         enabled = all(state_pn >= bmatrix' & ~any(inmatrix'>0 & inmatrix' <= state_pn));
                    %                     end

                    %prob_sync_p = sync{a}.passive{1}.prob(state_a, state_p)
                    %if prob_sync_p > 0
                    %% debugging block
                    %if options.verbose == 2
                    %    line_printf('---\n');
                    %    sync{a}.active{1}.print,
                    %    sync{a}.passive{1}.print
                    %end
                    %%
                    if node_p == node_a %self-loop
                        [new_state_p, ~, outprob_p] = State.afterEventHashed( sn, node_p, new_state_a(ia), event_p, class_p);
                    else % departure
                        [new_state_p, ~, outprob_p] = State.afterEventHashed( sn, node_p, state_p, event_p, class_p);
                    end
                    %  SPN code:
                    %                    if node_p == node_a %self-loop
                    %                         [new_state_p, ~, outprob_p, trans_p, modes_p] = State.afterEventHashed( qn, node_p, new_state_a(ia), event_p, class_p);
                    %                     else % departure
                    %                         [new_state_p, ~, outprob_p, trans_p, modes_p] = State.afterEventHashed( qn, node_p, state_p, event_p, class_p);
                    %                     end
                    for ip=1:size(new_state_p,1)
                        if node_p ~= local
                            if new_state_p ~= -1
                                if sn.isstatedep(node_a,3)
                                    newStateCell = stateCell;
                                    newStateCell{sn.nodeToStateful(node_a)} = sn.space{sn.nodeToStateful(node_a)}(new_state_a(ia),:);
                                    newStateCell{sn.nodeToStateful(node_p)} = sn.space{sn.nodeToStateful(node_p)}(new_state_p(ip),:);
                                    prob_sync_p = sync{a}.passive{1}.prob(stateCell, newStateCell) * outprob_p(ip); %state-dependent
                                else
                                    prob_sync_p = sync{a}.passive{1}.prob * outprob_p(ip);
                                end
                            else
                                prob_sync_p = 0;
                            end
                        end
                        if ~isempty(new_state_a(ia))
                            if node_p == local % local action
                                new_state = state;
                                new_state(sn.nodeToStateful(node_a)) = new_state_a(ia);
                                prob_sync_p = outprob_p(ip);
                            elseif ~isempty(new_state_p)
                                new_state = state;
                                new_state(sn.nodeToStateful(node_a)) = new_state_a(ia);
                                new_state(sn.nodeToStateful(node_p)) = new_state_p(ip);
                            end
                            % SPN code:
                            %       if enabled
                            %                                 ns = find(ismember(SSh(:,[sn.nodeToStateful(node_a),sn.nodeToStateful(node_p)]),[new_state_a(ia),new_state_p(ip)],'rows'));
                            %                                 for ins=1:length(ns)
                            %                                    if ns(ins) > 0 && ~isempty(trans_p)
                            %                                         tr = trans_p(ip);
                            %                                         mode = modes_p(ip);
                            %                                         bmatrix = sn.varsparam{tr}.back(:,mode);
                            %                                         fmatrix = sn.varsparam{tr}.forw(:,mode);
                            %                                         cmatrix = fmatrix - bmatrix;
                            %                                         if isequal(state_pn + cmatrix',SS(ns(ins),3:end))
                            %                                             [ex_a,seq_a] = ST.search(state_pn');
                            %                                             [ex_p,seq_p] = ST.search(SS(ns(ins),3:end)');
                            %                                             if ex_a && ex_p && edges(seq_a, seq_p)
                            %                                                 Adj_m(s, ns(ins)) = modes_p(ip);
                            %                                                 Adj_t(s, ns(ins)) = trans_p(ip);
                            % %                                                 s,ns(ins)
                            %                                                 if ~isnan(rate_a(ia))
                            %                                                     if node_p < local && ~csmask(class_a, class_p) && rate_a(ia) * prob_sync_p >0 && (sn.nodetype(node_p)~=NodeType.Source)
                            %                                                         error('Error: state-dependent routing at node %d (%s) violates the class switching mask (node %d -> node %d, class %d -> class %d).', node_a, sn.nodenames{node_a}, node_a, node_p, class_a, class_p);
                            %                                                     end
                            %                                                     if size(Dfilt{a}) >= [s,ns(ins)] % check needed as D{a} is a sparse matrix
                            %                                                         Dfilt{a}(s,ns(ins)) = Dfilt{a}(s,ns(ins)) + rate_a(ia) * prob_sync_p;
                            %                                                     else
                            %                                                         Dfilt{a}(s,ns(ins)) = rate_a(ia) * prob_sync_p;
                            %                                                     end
                            %                                                 end
                            %                                             end
                            %                                         end
                            %                                    end
                            %                                 end
                            %            else

                            ns = matchrow(stateSpaceHashed, new_state);
                            if ns>0
                                if ~isnan(rate_a)
                                    if node_p < local && ~csmask(class_a, class_p) && rate_a(ia) * prob_sync_p >0 && (sn.nodetype(node_p)~=NodeType.Source)
                                        line_error(mfilename,sprintf('Error: state-dependent routing at node %d (%s) violates the class switching mask (node %s -> node %s, class %s -> class %s).', node_a, sn.nodenames{node_a}, sn.nodenames{node_a}, sn.nodenames{node_p}, sn.classnames{class_a}, sn.classnames{class_p}));
                                    end
                                    if size(Dfilt{a}) >= [s,ns] % check needed as D{a} is a sparse matrix
                                        Dfilt{a}(s,ns) = Dfilt{a}(s,ns) + rate_a(ia) * prob_sync_p;
                                    else
                                        Dfilt{a}(s,ns) = rate_a(ia) * prob_sync_p;
                                    end
                                end
                            end
                            % SPN code:
                            %                            end
                        end
                    end
                else % node_p == local
                    if ~isempty(new_state_a(ia))
                        new_state = state;
                        new_state(sn.nodeToStateful(node_a)) = new_state_a(ia);
                        prob_sync_p = 1;
                        ns = matchrow(stateSpaceHashed, new_state);
                        if ns>0
                            if ~isnan(rate_a)
                                if size(Dfilt{a}) >= [s,ns] % needed for sparse matrix
                                    Dfilt{a}(s,ns) = Dfilt{a}(s,ns) + rate_a(ia) * prob_sync_p;
                                else
                                    Dfilt{a}(s,ns) = rate_a(ia) * prob_sync_p;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

for a=1:A
    Q = Q + Dfilt{a};
end
Q = Q - diag(diag(Q));
%SolverCTMC.printInfGen(Q,stateSpace)
%%
arvRates = zeros(size(stateSpaceHashed,1),nstateful,nclasses);
depRates = zeros(size(stateSpaceHashed,1),nstateful,nclasses);
for a=1:A
    % active
    node_a = sync{a}.active{1}.node;
    class_a = sync{a}.active{1}.class;
    event_a = sync{a}.active{1}.event;
    % passive
    node_p = sync{a}.passive{1}.node;
    class_p = sync{a}.passive{1}.class;
    if event_a == EventType.DEP
        node_a_sf = sn.nodeToStateful(node_a);
        node_p_sf = sn.nodeToStateful(node_p);
        for s=1:size(stateSpaceHashed,1)
            depRates(s,node_a_sf,class_a) = depRates(s,node_a_sf,class_a) + sum(Dfilt{a}(s,:));
            arvRates(s,node_p_sf,class_p) = arvRates(s,node_p_sf,class_p) + sum(Dfilt{a}(s,:));
        end
    end
end
zero_row = find(sum(Q,2)==0);
zero_col = find(sum(Q,1)==0);

%%
% in case the last column of Q represent a state for a transient class, it
% is possible that no transitions go back to it, although it is valid for
% the system to be initialized in that state. So we need to fill-in the
% zeros at the end.
Q(:,end+1:end+(size(Q,1)-size(Q,2)))=0;
Q(zero_row,zero_row) = -eye(length(zero_row)); % can this be replaced by []?
Q(zero_col,zero_col) = -eye(length(zero_col));
for a=1:A
    Dfilt{a}(:,end+1:end+(size(Dfilt{a},1)-size(Dfilt{a},2)))=0;
end

if options.verbose == VerboseLevel.DEBUG || GlobalConstants.Verbose == VerboseLevel.DEBUG
    SolverCTMC.printInfGen(Q,stateSpace);
end
Q = ctmc_makeinfgen(Q);
%% now remove immediate transitions
% we first determine states in stateful nodes where there is an immediate
% job in the node

if options.config.hide_immediate % if want to remove immediate transitions
    imm = [];
    for ind = 1:sn.nnodes
        % Only consider non-station stateful nodes that are not Caches
        % Cache nodes need their immediate transitions to compute hit/miss rates
        if sn.isstateful(ind) && ~sn.isstation(ind) && sn.nodetype(ind) ~= NodeType.Cache
            isf = sn.nodeToStateful(ind);
            imm_st = find(sum(sn.space{isf}(:,1:nclasses),2)>0);
            imm = [imm; find(arrayfun(@(a) any(a==imm_st),stateSpaceHashed(:,isf)))];
        end
    end
    imm = unique(imm);
    nonimm = setdiff(1:size(Q,1),imm);
    stateSpace(imm,:) = [];
    stateSpaceAggr(imm,:) = [];
    %    full(Q)
    [Q,~,Q12,~,Q22] = ctmc_stochcomp(Q, nonimm);
    %    full(Q)
    for a=1:A
        % stochastic complement for action a
        Q21a = Dfilt{a}(imm,nonimm);
        Ta = (-Q22) \ Q21a;
        Ta = Q12*Ta;
        Dfilt{a} = Dfilt{a}(nonimm,nonimm)+Ta;
    end

    depRates(imm,:,:) = [];
    arvRates(imm,:,:) = [];
    % recompute arvRates and depRates
    %     arvRates = zeros(size(stateSpace,1),nstateful,nclasses);
    %     depRates = zeros(size(stateSpace,1),nstateful,nclasses);
    %     for a=1:A
    %         % active
    %         node_a = sync{a}.active{1}.node;
    %         class_a = sync{a}.active{1}.class;
    %         event_a = sync{a}.active{1}.event;
    %         % passive
    %         node_p = sync{a}.passive{1}.node;
    %         class_p = sync{a}.passive{1}.class;
    %         if event_a == EventType.DEP
    %             node_a_sf = sn.nodeToStateful(node_a);
    %             node_p_sf = sn.nodeToStateful(node_p);
    %             for s=1:size(stateSpace,1)
    %                 depRates(s,node_a_sf,class_a) = depRates(s,node_a_sf,class_a) + sum(Dfilt{a}(s,:));
    %                 arvRates(s,node_p_sf,class_p) = arvRates(s,node_p_sf,class_p) + sum(Dfilt{a}(s,:));
    %             end
    %         end
    %     end
end
%SolverCTMC.printInfGen(Q,stateSpace)
%
% Draft SPN:
% if ~isempty(Adj) && ~isempty(ST)
%     imm = [];
%     for s=1:size(SSh, 1)
%         % check for immediate transitions
%         enabled_t = find(Adj_t(s,:));
%         imm_t = [];
%         for t=1:length(enabled_t)
%             if sn.varsparam{Adj_t(s,enabled_t(t))}.timingstrategies(Adj_m(s,enabled_t(t)))
%                 imm_t = [imm_t, enabled_t(t)];
%             end
%         end
%         % Immediate transitions exists, the marking needs to be removed.
%         if ~isempty(imm_t)
%             imm = [imm, s];
%             non_imm = setdiff(enabled_t,imm_t);
%             inM = find(Adj_t(:,s));
%             for in=1:length(inM)
%                 depRates(inM(in),:,:) = depRates(inM(in),:,:) + depRates(s,:,:);
%                 for nim=1:length(non_imm)
%                     Q(inM(in), non_imm(nim)) = Q(inM(in), non_imm(nim)) + Q(inM(in), s) + Q(s, non_imm(nim));
%                     if Adj_t(inM(in), non_imm(nim))==0
%                         Adj_t(inM(in), non_imm(nim)) = Adj_t(inM(in), s);
%                         Adj_m(inM(in), non_imm(nim)) = Adj_m(inM(in), s);
%                     end
%                     arvRates(non_imm(nim),:,:) = arvRates(non_imm(nim),:,:) + arvRates(s,:,:);
%                 end
%                 totalWeight = 0;
%                 for im=1:length(imm_t)
%                     weight = sn.varsparam{Adj_t(s,imm_t(im))}.firingweights(Adj_m(s,imm_t(im)));
%                     totalWeight = totalWeight + weight;
%                     if Adj_t(inM(in), imm_t(im))==0
%                         Adj_t(inM(in), imm_t(im)) = Adj_t(inM(in), s);
%                         Adj_m(inM(in), imm_t(im)) = Adj_m(inM(in), s);
%                     end
%                 end
%                 for im=1:length(imm_t)
%                     weight = sn.varsparam{Adj_t(s,imm_t(im))}.firingweights(Adj_m(s,imm_t(im)));
%                     Q(inM(in), imm_t(im)) = Q(inM(in), imm_t(im)) + Q(inM(in), s) * (weight/totalWeight);
%                 end
%             end
%             Adj_t(s,:) = 0;
%             Adj_m(s,:) = 0;
%             Adj_t(:,s) = 0;
%             Adj_m(:,s) = 0;
%         end
%     end
%     Q(imm,:) = [];
%     SS(imm,:) = [];
%     SSq(imm,:) = [];
%     arvRates(imm,:,:) = [];
%     depRates(imm,:,:) = [];
%     Q(:,imm) = [];
% end
%%
%Q = ctmc_makeinfgen(Q);
end
