function [nodeStateSpace, sn, capacityc] = spaceGeneratorNodes(sn, cutoff, options)
if nargin<3
    options = Solver.defaultOptions;
end
N = sn.njobs';
sn.space = {};
capacityc = zeros(sn.nnodes, sn.nclasses);
% Draft SPN support
% for n=1:sn.nnodes
%     if isfield(sn.varsparam{n}, "capacityc")
%         capacityc(n) = sn.varsparam{n}.capacityc;
%     end
% end
for ind=1:sn.nnodes
    if sn.isstation(ind) % place jobs across stations
        ist = sn.nodeToStation(ind);
        isf = sn.nodeToStateful(ind);
        for r=1:sn.nclasses %cut-off open classes to finite capacity
            c = find(sn.chains(:,r));
            if ~isempty(sn.visits{c}) && sn.visits{c}(isf,r) == 0
                capacityc(ind,r) = 0;
            % Draft SPN support
            %elseif isfield(sn.varsparam{ind}, 'capacityc') 
            %    capacityc(ind,r) =  min(cutoff(ist,r), sn.classcap(ist,r));                
            elseif ~isempty(sn.proc) && ~isempty(sn.proc{ist}{r}) && any(any(isnan(sn.proc{ist}{r}{1}))) % disabled
                capacityc(ind,r) = 0;
            else
                if isinf(N(r))
                    capacityc(ind,r) =  min(cutoff(ist,r), sn.classcap(ist,r));
                else
                    capacityc(ind,r) =  sum(sn.njobs(sn.chains(c,:)));
                end
            end
        end
        if sn.isstation(ind)
            % in this case, the local variables are produced within
            % fromMarginalBounds, e.g., for RROBIN routing
            sn.space{isf} = State.fromMarginalBounds(sn, ind, [], capacityc(ind,:), sn.cap(ist), options);
        else
            % this is the case for example of cache nodes
            state_bufsrv = State.fromMarginalBounds(sn, ind, [], capacityc(ind,:), sn.cap(ist), options);
            state_var = State.spaceLocalVars(sn, ind);
            sn.space{isf} = State.cartesian(state_bufsrv,state_var); % generate all possible states for local variables
        end
        if isinf(sn.nservers(ist))
            sn.nservers(ist) = sum(capacityc(ind,:));
        end
    elseif sn.isstateful(ind) % generate state space of other stateful nodes that are not stations
        %ist = sn.nodeToStation(ind);
        isf = sn.nodeToStateful(ind);
        switch sn.nodetype(ind)
            case NodeType.Cache
                for r=1:sn.nclasses % restrict state space generation to immediate events
                    if isnan(sn.nodeparam{ind}.pread{r})
                        capacityc(ind,r) =  1; %
                    else
                        capacityc(ind,r) =  1; %
                    end
                end
            case NodeType.Router
                % For Router nodes, only allow capacity for classes that have
                % non-zero nodevisits (classes that actually visit the Router)
                for r=1:sn.nclasses
                    c = find(sn.chains(:,r));
                    if ~isempty(sn.nodevisits{c}) && sn.nodevisits{c}(ind,r) > 0
                        capacityc(ind,r) = 1;
                    else
                        capacityc(ind,r) = 0;
                    end
                end
            otherwise
                capacityc(ind,:) =  1; %
        end
        state_bufsrv = State.fromMarginalBounds(sn, ind, [], capacityc(ind,:), 1, options);
        state_var = State.spaceLocalVars(sn, ind);
        sn.space{isf} = State.cartesian(state_bufsrv,state_var); % generate all possible states for local variables
    end
end
nodeStateSpace = sn.space;
end