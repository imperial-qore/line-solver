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
    % Cooperative wall-clock budget checkpoint (options.timeout)
    if nargin>=3 && lineTimeoutExceeded(options)
        line_error(mfilename,'State space generation exceeded the wall-clock time budget (options.timeout=%gs).', options.timeout);
    end
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
            elseif sn.nodetype(ind) ~= NodeType.Place && ~isempty(sn.proc) && ~isempty(sn.proc{ist}{r}) && any(any(isnan(sn.proc{ist}{r}{1}))) % disabled (not Places - they hold tokens without service)
                capacityc(ind,r) = 0;
            else
                if isinf(N(r))
                    capacityc(ind,r) =  min(cutoff(ist,r), sn.classcap(ist,r));
                else
                    % closed classes: enumerate up to the chain population, but never
                    % beyond the class capacity at this station (finite-buffer stations)
                    capacityc(ind,r) =  min(sum(sn.njobs(sn.chains(c,:))), sn.classcap(ist,r));
                end
            end
            % Finite-capacity region bound: a station in a DROP/WAITQ region can
            % never hold more than the region's per-class cap (nor more than the
            % region-global cap) of class r, so enumerating beyond it produces only
            % states the region filter later discards. Bounding capacityc here keeps
            % the generated space small (the auto-cutoff is region-blind and can be
            % far larger than any reachable population), with an identical final
            % result. -1 = unbounded. sn.region{f} is M x (K+1): cols 1..K per-class,
            % col K+1 the region-global cap at each member station.
            if isfield(sn,'nregions') && sn.nregions > 0 && capacityc(ind,r) > 0
                for f=1:sn.nregions
                    regf = sn.region{f};
                    if ist <= size(regf,1)
                        rc = regf(ist, r);
                        gc = regf(ist, sn.nclasses+1);
                        if rc >= 0
                            capacityc(ind,r) = min(capacityc(ind,r), rc);
                        end
                        if gc >= 0
                            capacityc(ind,r) = min(capacityc(ind,r), gc);
                        end
                    end
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
                    capacityc(ind,r) = 1;
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
            case NodeType.Transition
                capacityc(ind,:) = 0; % Transitions don't hold class-based jobs
                % Generate per-mode state space (bypass fromMarginalBounds)
                nmodes = sn.nodeparam{ind}.nmodes;
                firingphases = sn.nodeparam{ind}.firingphases;
                if any(isnan(firingphases))
                    firingphases = zeros(1, nmodes);
                    for m = 1:nmodes
                        if iscell(sn.nodeparam{ind}.firingproc) && ~isempty(sn.nodeparam{ind}.firingproc{m})
                            firingphases(m) = size(sn.nodeparam{ind}.firingproc{m}{1}, 1);
                        else
                            firingphases(m) = 1;
                        end
                    end
                end
                fK = firingphases;
                nmodeservers = sn.nodeparam{ind}.nmodeservers;
                max_jobs = sum(sn.njobs(~isinf(sn.njobs)));
                if any(isinf(sn.njobs))
                    max_jobs = max_jobs + sum(cutoff(1,:));
                end
                mode_spaces = cell(1, nmodes);
                for m = 1:nmodes
                    max_srv_m = nmodeservers(m);
                    if isinf(max_srv_m)
                        max_srv_m = max_jobs;
                    end
                    max_srv_m = min(max_srv_m, max_jobs);
                    mode_states = [];
                    for total = 0:max_srv_m
                        phase_combs = multichoose(fK(m), total);
                        buf_m = nmodeservers(m);
                        if isinf(buf_m)
                            buf_m = GlobalConstants.MaxInt();
                        end
                        buf_m = buf_m - total;
                        mode_states = [mode_states; repmat(buf_m, size(phase_combs,1), 1), phase_combs]; %#ok<AGROW>
                    end
                    mode_spaces{m} = mode_states;
                end
                trans_space = mode_spaces{1};
                for m = 2:nmodes
                    trans_space = State.cartesian(trans_space, mode_spaces{m});
                end
                % State.cartesian interleaves the columns as
                % [idle_1, phases_1, idle_2, phases_2, ...]. Reorder to the
                % contiguous layout [idle(nmodes), phases(sum fK), fired(nmodes)]
                % expected by afterGlobalEvent.m (space_buf = 1:nmodes,
                % space_srv = nmodes+(1:sum fK)) and produced by the native-Python
                % builder. The interleaved layout made afterGlobalEvent read the
                % idle/phase counts of every mode after the first from the wrong
                % columns, which collapsed the state space of any transition with
                % more than one mode.
                idle_cols = zeros(1, nmodes);
                phase_cols = [];
                off = 0;
                for m = 1:nmodes
                    idle_cols(m) = off + 1;
                    phase_cols = [phase_cols, (off+2):(off+1+fK(m))]; %#ok<AGROW>
                    off = off + 1 + fK(m);
                end
                trans_space = trans_space(:, [idle_cols, phase_cols]);
                trans_space = [trans_space, zeros(size(trans_space,1), nmodes)]; % append fired counts
                state_var = State.spaceLocalVars(sn, ind);
                sn.space{isf} = State.cartesian(trans_space, state_var);
                continue; % skip fromMarginalBounds below
            otherwise
                capacityc(ind,:) =  1; %
        end
        state_bufsrv = State.fromMarginalBounds(sn, ind, [], capacityc(ind,:), 1, options);
        state_var = State.spaceLocalVars(sn, ind);
        sn.space{isf} = State.cartesian(state_bufsrv,state_var); % generate all possible states for local variables

        % Prune locally invalid cache states for retrieval-system caches. The
        % unfiltered cartesian product enumerates states the dynamics can never
        % reach -- an item simultaneously cached and being retrieved, more items
        % in retrieval than the retrieval-system capacity, or more than one job
        % in the cache server. Leaving them in pollutes the CTMC stationary
        % distribution. The cross-node (cache <-> queue) consistency is enforced
        % during global composition. Non-retrieval caches are left untouched.
        if sn.nodetype(ind) == NodeType.Cache ...
                && isfield(sn.nodeparam{ind},'retrievalSystemCapacity') ...
                && sn.nodeparam{ind}.retrievalSystemCapacity > 0
            np = sn.nodeparam{ind};
            nItems = np.nitems;
            tcc = np.totalCacheCapacity;
            rsc = np.retrievalSystemCapacity;
            lvs = size(state_bufsrv,2);          % per-class server presence is columns 1..lvs
            value = sn.space{isf};
            keep = false(size(value,1),1);
            for row = 1:size(value,1)
                st = value(row,:);
                % at most one class job may sit in the cache server/buffer
                validRow = sum(st(1:lvs)) <= 1;
                if validRow
                    nbits = 0;
                    cacheItems = st((lvs+1):(lvs+tcc));
                    for col = (lvs+tcc+1):size(st,2)
                        if st(col) == 0, continue; end
                        nbits = nbits + 1;
                        item = col - (lvs+tcc);  % 1-based item id of the set bit
                        % a cached item cannot simultaneously be in the retrieval system
                        if any(cacheItems == item)
                            validRow = false; break;
                        end
                    end
                    % at most retrievalSystemCapacity items may be in retrieval at once
                    if validRow && nbits > rsc
                        validRow = false;
                    end
                end
                keep(row) = validRow;
            end
            sn.space{isf} = value(keep,:);
        end
    end
end
nodeStateSpace = sn.space;
end