classdef JLINE
    % JLINE Conversion utilities for JLINE format models
    %
    % JLINE provides static methods to convert between LINE MATLAB models 
    % and JLINE Java/Kotlin models. This class serves as the primary interface
    % for interoperability between the MATLAB and Java implementations of LINE.
    %
    % @brief JLINE format conversion and Java interoperability utilities
    %
    % Main functionality:
    % - Convert LINE MATLAB models to JLINE Java models
    % - Convert JLINE Java models back to LINE MATLAB format  
    % - Access JLINE solvers from MATLAB
    % - Handle serialization between MATLAB and Java representations
    %
    % Example:
    % @code
    % % Convert a LINE model to JLINE format
    % jnetwork = JLINE.from_model(network);
    % % Get a JLINE solver
    % jssa = JLINE.get_solver(jnetwork, 'ssa');
    % @endcode

    properties (Constant)
        jar_loc = which('jline.jar');
    end

    methods(Static)

        function model = from_line_layered_network(line_layered_network)
            sn = line_layered_network.getStruct;

            %% initialization
            model = jline.lang.layered.LayeredNetwork(line_layered_network.getName);

            %% host processors
            P = cell(1,sn.nhosts);
            for h=1:sn.nhosts
                if isinf(sn.mult(h))
                    sn_mult_h = java.lang.Integer.MAX_VALUE;
                else
                    sn_mult_h = sn.mult(h);
                end
                switch sn.sched(h)
                    case SchedStrategy.REF
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.REF);
                    case SchedStrategy.INF
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.INF);
                    case SchedStrategy.FCFS
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.FCFS);
                    case SchedStrategy.LCFS
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFS);
                    case SchedStrategy.SIRO
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.SIRO);
                    case SchedStrategy.SJF
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.SJF);
                    case SchedStrategy.LJF
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LJF);
                    case SchedStrategy.PS
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.PS);
                    case SchedStrategy.DPS
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.DPS);
                    case SchedStrategy.GPS
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.GPS);
                    case SchedStrategy.SEPT
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.SEPT);
                    case SchedStrategy.LEPT
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LEPT);
                    case {SchedStrategy.HOL, SchedStrategy.FCFSPRIO}
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.FCFSPRIO);
                    case SchedStrategy.FORK
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.FORK);
                    case SchedStrategy.EXT
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.EXT);
                    case SchedStrategy.LCFSPR
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFSPR);
                    case SchedStrategy.LCFSPI
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFSPI);
                    case SchedStrategy.LCFSPRIO
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFSPRIO);
                    case SchedStrategy.LCFSPRPRIO
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFSPRPRIO);
                    case SchedStrategy.LCFSPIPRIO
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFSPIPRIO);
                    case SchedStrategy.PSPRIO % todo
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.PSPRIO);
                    case SchedStrategy.DPSPRIO % todo
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.DPSPRIO);
                    case SchedStrategy.GPSPRIO % todo
                        P{h} = jline.lang.layered.Processor(model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.GPSPRIO);
                end
                if sn.repl(h)~=1
                    P{h}.setReplication(sn.repl(h));
                end
            end

            %% tasks
            T = cell(1,sn.ntasks);
            for t=1:sn.ntasks
                tidx = sn.tshift+t;
                if isinf(sn.mult(tidx))
                    sn_mult_tidx = java.lang.Integer.MAX_VALUE;
                else
                    sn_mult_tidx = sn.mult(tidx);
                end
                % Check if this is a CacheTask
                if sn.iscache(tidx)
                    % Get replacement strategy from ordinal
                    switch sn.replacestrat(tidx)
                        case ReplacementStrategy.RR
                            jReplacestrat = jline.lang.constant.ReplacementStrategy.RR;
                        case ReplacementStrategy.FIFO
                            jReplacestrat = jline.lang.constant.ReplacementStrategy.FIFO;
                        case ReplacementStrategy.SFIFO
                            jReplacestrat = jline.lang.constant.ReplacementStrategy.SFIFO;
                        case ReplacementStrategy.LRU
                            jReplacestrat = jline.lang.constant.ReplacementStrategy.LRU;
                        otherwise
                            jReplacestrat = jline.lang.constant.ReplacementStrategy.FIFO;
                    end
                    T{t} = jline.lang.layered.CacheTask(model, sn.names{tidx}, sn.nitems(tidx), sn.itemcap{tidx}, jReplacestrat, sn_mult_tidx);
                else
                    switch sn.sched(tidx)
                        case SchedStrategy.REF
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.REF);
                        case SchedStrategy.INF
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.INF);
                        case SchedStrategy.FCFS
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.FCFS);
                        case SchedStrategy.LCFS
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.LCFS);
                        case SchedStrategy.SIRO
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.SIRO);
                        case SchedStrategy.SJF
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.SJF);
                        case SchedStrategy.LJF
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.LJF);
                        case SchedStrategy.PS
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.PS);
                        case SchedStrategy.DPS
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.DPS);
                        case SchedStrategy.GPS
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.GPS);
                        case SchedStrategy.SEPT
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.SEPT);
                        case SchedStrategy.LEPT
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.LEPT);
                        case {SchedStrategy.HOL, SchedStrategy.FCFSPRIO}
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.FCFSPRIO);
                        case SchedStrategy.FORK
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.FORK);
                        case SchedStrategy.EXT
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.EXT);
                        case SchedStrategy.LCFSPR
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.LCFSPR);
                        case SchedStrategy.PSPRIO % todo
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.PSPRIO);
                        case SchedStrategy.DPSPRIO % todo
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.DPSPRIO);
                        case SchedStrategy.GPSPRIO % todo
                            T{t} = jline.lang.layered.Task(model, sn.names{tidx}, sn_mult_tidx, jline.lang.constant.SchedStrategy.GPSPRIO);
                    end
                end
                T{t}.on(P{sn.parent(tidx)});
                if sn.repl(tidx)~=1
                    T{t}.setReplication(sn.repl(tidx));
                end
                if ~isempty(sn.think{tidx}) && sn.think_type(tidx) ~= ProcessType.DISABLED
                    switch sn.think_type(tidx)
                        case ProcessType.IMMEDIATE
                            T{t}.setThinkTime(jline.lang.processes.Immediate);
                        case ProcessType.EXP
                            T{t}.setThinkTime(jline.lang.processes.Exp(1/sn.think_mean(tidx)));
                        case ProcessType.ERLANG
                            T{t}.setThinkTime(jline.lang.processes.Erlang.fitMeanAndSCV(sn.think_mean(tidx), sn.think_scv(tidx)));
                        case ProcessType.HYPEREXP
                            if ~isempty(sn.think_params{tidx}) && length(sn.think_params{tidx}) >= 3
                                T{t}.setThinkTime(jline.lang.processes.HyperExp(sn.think_params{tidx}(1), sn.think_params{tidx}(2), sn.think_params{tidx}(3)));
                            else
                                T{t}.setThinkTime(jline.lang.processes.HyperExp.fitMeanAndSCV(sn.think_mean(tidx), sn.think_scv(tidx)));
                            end
                        case ProcessType.COXIAN
                            T{t}.setThinkTime(jline.lang.processes.Coxian.fitMeanAndSCV(sn.think_mean(tidx), sn.think_scv(tidx)));
                        case ProcessType.APH
                            T{t}.setThinkTime(jline.lang.processes.APH.fitMeanAndSCV(sn.think_mean(tidx), sn.think_scv(tidx)));
                        case ProcessType.PH
                            if ~isempty(sn.think_proc{tidx})
                                proc = sn.think_proc{tidx};
                                T{t}.setThinkTime(jline.lang.processes.PH(JLINE.from_line_matrix(proc{1}), JLINE.from_line_matrix(proc{2})));
                            else
                                T{t}.setThinkTime(jline.lang.processes.Exp(1/sn.think_mean(tidx)));
                            end
                        case ProcessType.MAP
                            if ~isempty(sn.think_proc{tidx})
                                proc = sn.think_proc{tidx};
                                T{t}.setThinkTime(jline.lang.processes.MAP(JLINE.from_line_matrix(proc{1}), JLINE.from_line_matrix(proc{2})));
                            else
                                T{t}.setThinkTime(jline.lang.processes.Exp(1/sn.think_mean(tidx)));
                            end
                        otherwise
                            line_error(mfilename,sprintf('JLINE conversion does not support the %s distribution for task think time yet.',char(sn.think_type(tidx))));
                    end
                end
            end
            %% entries
            E = cell(1,sn.nentries);
            for e=1:sn.nentries
                eidx = sn.eshift+e;
                % Check if this is an ItemEntry (has nitems > 0)
                if sn.nitems(eidx) > 0
                    % ItemEntry requires cardinality and popularity distribution
                    if ~isempty(sn.itemproc) && ~isempty(sn.itemproc{eidx})
                        jPopularity = JLINE.from_line_distribution(sn.itemproc{eidx});
                    else
                        % Default to uniform distribution
                        jPopularity = jline.lang.processes.DiscreteSampler(jline.util.matrix.Matrix.uniformDistribution(sn.nitems(eidx)));
                    end
                    E{e} = jline.lang.layered.ItemEntry(model, sn.names{eidx}, sn.nitems(eidx), jPopularity);
                else
                    E{e} = jline.lang.layered.Entry(model, sn.names{eidx});
                end
                E{e}.on(T{sn.parent(eidx)-sn.tshift});
            end

            %% activities
            A = cell(1,sn.nacts);
            for a=1:sn.nacts
                aidx = sn.ashift+a;
                tidx = sn.parent(aidx);
                onTask = tidx-sn.tshift;
                % Convert host demand from primitives to Java distribution
                switch sn.hostdem_type(aidx)
                    case ProcessType.IMMEDIATE
                        jHostDem = jline.lang.processes.Immediate;
                    case ProcessType.DISABLED
                        jHostDem = jline.lang.processes.Disabled;
                    case ProcessType.EXP
                        jHostDem = jline.lang.processes.Exp(1/sn.hostdem_mean(aidx));
                    case ProcessType.ERLANG
                        jHostDem = jline.lang.processes.Erlang.fitMeanAndSCV(sn.hostdem_mean(aidx), sn.hostdem_scv(aidx));
                    case ProcessType.HYPEREXP
                        if ~isempty(sn.hostdem_params{aidx}) && length(sn.hostdem_params{aidx}) >= 3
                            jHostDem = jline.lang.processes.HyperExp(sn.hostdem_params{aidx}(1), sn.hostdem_params{aidx}(2), sn.hostdem_params{aidx}(3));
                        else
                            jHostDem = jline.lang.processes.HyperExp.fitMeanAndSCV(sn.hostdem_mean(aidx), sn.hostdem_scv(aidx));
                        end
                    case ProcessType.COXIAN
                        jHostDem = jline.lang.processes.Coxian.fitMeanAndSCV(sn.hostdem_mean(aidx), sn.hostdem_scv(aidx));
                    case ProcessType.APH
                        jHostDem = jline.lang.processes.APH.fitMeanAndSCV(sn.hostdem_mean(aidx), sn.hostdem_scv(aidx));
                    case ProcessType.PH
                        if ~isempty(sn.hostdem_proc{aidx})
                            proc = sn.hostdem_proc{aidx};
                            jHostDem = jline.lang.processes.PH(JLINE.from_line_matrix(proc{1}), JLINE.from_line_matrix(proc{2}));
                        else
                            jHostDem = jline.lang.processes.Exp(1/sn.hostdem_mean(aidx));
                        end
                    case ProcessType.MAP
                        if ~isempty(sn.hostdem_proc{aidx})
                            proc = sn.hostdem_proc{aidx};
                            jHostDem = jline.lang.processes.MAP(JLINE.from_line_matrix(proc{1}), JLINE.from_line_matrix(proc{2}));
                        else
                            jHostDem = jline.lang.processes.Exp(1/sn.hostdem_mean(aidx));
                        end
                    case ProcessType.DET
                        jHostDem = jline.lang.processes.Det(sn.hostdem_mean(aidx));
                    otherwise
                        line_error(mfilename,sprintf('JLINE conversion does not support the %s distribution for host demand yet.',char(sn.hostdem_type(aidx))));
                end
                A{a} = jline.lang.layered.Activity(model, sn.names{aidx}, jHostDem);
                A{a}.on(T{onTask});

                boundTo = find(sn.graph((sn.eshift+1):(sn.eshift+sn.nentries),aidx));

                if ~isempty(boundTo)
                    A{a}.boundTo(E{boundTo});
                end

                if sn.sched(tidx) ~= SchedStrategy.REF % ref tasks don't reply
                    repliesTo = find(sn.replygraph(a,:)); % index of entry
                    if ~isempty(repliesTo)
                        if ~sn.isref(sn.parent(sn.eshift+repliesTo))
                            A{a}.repliesTo(E{repliesTo});
                        end
                    end
                end

                if ~isempty(sn.callpair)
                    cidxs = find(sn.callpair(:,1)==aidx);
                    calls = sn.callpair(:,2);
                    for c = cidxs(:)'
                        switch sn.calltype(c)
                            case CallType.SYNC
                                A{a}.synchCall(E{calls(c)-sn.eshift},sn.callproc_mean(c));
                            case CallType.ASYNC
                                A{a}.asynchCall(E{calls(c)-sn.eshift},sn.callproc_mean(c));
                        end
                    end
                end

            end

            %% think times
            for h=1:sn.nhosts
                if ~isempty(sn.think{h}) && sn.think_type(h) ~= ProcessType.DISABLED
                    switch sn.think_type(h)
                        case ProcessType.IMMEDIATE
                            P{h}.setThinkTime(jline.lang.processes.Immediate);
                        case ProcessType.EXP
                            P{h}.setThinkTime(jline.lang.processes.Exp(1/sn.think_mean(h)));
                        case ProcessType.ERLANG
                            P{h}.setThinkTime(jline.lang.processes.Erlang.fitMeanAndSCV(sn.think_mean(h),sn.think_scv(h)));
                        case ProcessType.HYPEREXP
                            % For HyperExp, reconstruct from params if available, otherwise use fitMeanAndSCV
                            if ~isempty(sn.think_params{h}) && length(sn.think_params{h}) >= 3
                                P{h}.setThinkTime(jline.lang.processes.HyperExp(sn.think_params{h}(1), sn.think_params{h}(2), sn.think_params{h}(3)));
                            else
                                P{h}.setThinkTime(jline.lang.processes.HyperExp.fitMeanAndSCV(sn.think_mean(h), sn.think_scv(h)));
                            end
                        case ProcessType.COXIAN
                            % For Coxian, use fitMeanAndSCV
                            P{h}.setThinkTime(jline.lang.processes.Coxian.fitMeanAndSCV(sn.think_mean(h), sn.think_scv(h)));
                        case ProcessType.APH
                            % For APH, reconstruct from params if available
                            if ~isempty(sn.think_params{h})
                                P{h}.setThinkTime(jline.lang.processes.APH.fitMeanAndSCV(sn.think_mean(h), sn.think_scv(h)));
                            else
                                P{h}.setThinkTime(jline.lang.processes.Exp(1/sn.think_mean(h)));
                            end
                        otherwise
                            line_error(mfilename,sprintf('JLINE conversion does not support the %s distribution yet.',char(sn.think_type(h))));
                    end
                end
            end

            %% Sequential precedences
            for ai = 1:sn.nacts
                aidx = sn.ashift + ai;
                tidx = sn.parent(aidx);
                % for all successors
                for bidx=find(sn.graph(aidx,:))
                    if bidx > sn.ashift % ignore precedence between entries and activities
                        % Serial pattern (SEQ)
                        if full(sn.actpretype(aidx)) == ActivityPrecedenceType.PRE_SEQ && full(sn.actposttype(bidx)) == ActivityPrecedenceType.POST_SEQ
                            T{tidx-sn.tshift}.addPrecedence(jline.lang.layered.ActivityPrecedence.Serial(sn.names{aidx}, sn.names{bidx}));
                        end
                    end
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%% translated up to here

            %% Loop precedences (POST_LOOP)
            % Loop structure in sn.graph:
            %   - Entry activity (preAct) has edge to first loop body activity
            %   - Last loop body activity has back-edge to first loop body (weight = 1-1/counts)
            %   - Last loop body activity has edge to end activity (weight = 1/counts)
            %   - All loop body and end activities have POST_LOOP type
            processedLoops = false(1, sn.nacts);
            for ai = 1:sn.nacts
                aidx = sn.ashift + ai;
                tidx = sn.parent(aidx);
                % Check if this activity starts a loop (has a successor with POST_LOOP type)
                % and hasn't been processed as part of another loop
                if processedLoops(ai)
                    continue;
                end

                successors = find(sn.graph(aidx,:));
                for bidx = successors
                    if bidx > sn.ashift && full(sn.actposttype(bidx)) == ActivityPrecedenceType.POST_LOOP
                        % Skip if this loop body activity was already processed
                        if processedLoops(bidx - sn.ashift)
                            continue;
                        end
                        % Found start of a loop: aidx is the entry, bidx is first loop body activity
                        loopStart = bidx;
                        precActs = java.util.ArrayList();

                        % Follow the chain of POST_LOOP activities
                        curIdx = loopStart;
                        while true
                            precActs.add(sprintf("%s", sn.names{curIdx}));
                            processedLoops(curIdx - sn.ashift) = true;

                            % Find successors of current activity
                            curSuccessors = find(sn.graph(curIdx,:));
                            curSuccessors = curSuccessors(curSuccessors > sn.ashift);

                            % Check for loop termination: find the end activity
                            % End activity has weight = 1/counts (not the back-edge weight)
                            endIdx = 0;
                            nextIdx = 0;
                            for succIdx = curSuccessors
                                if full(sn.actposttype(succIdx)) == ActivityPrecedenceType.POST_LOOP
                                    if succIdx == loopStart
                                        % This is the back-edge, skip it
                                        continue;
                                    end
                                    weight = full(sn.graph(curIdx, succIdx));
                                    if weight > 0 && weight < 1
                                        % This is the end activity (weight = 1/counts)
                                        endIdx = succIdx;
                                    else
                                        % This is the next activity in the loop body (weight = 1.0)
                                        nextIdx = succIdx;
                                    end
                                end
                            end

                            if endIdx > 0
                                % Found end activity - calculate counts and output
                                weight = full(sn.graph(curIdx, endIdx));
                                if weight > 0
                                    counts = 1/weight;
                                else
                                    counts = 1; % Fallback to prevent division by zero
                                end
                                precActs.add(sprintf("%s", sn.names{endIdx}));
                                processedLoops(endIdx - sn.ashift) = true;

                                T{tidx-sn.tshift}.addPrecedence(jline.lang.layered.ActivityPrecedence.Loop(sn.names{aidx}, precActs, jline.util.matrix.Matrix(counts)));
                                break;
                            elseif nextIdx > 0
                                % Continue to next activity in loop body
                                curIdx = nextIdx;
                            else
                                % No more successors - shouldn't happen in valid loop
                                break;
                            end
                        end
                        break; % Only process one loop starting from this activity
                    end
                end
            end

            %% OrFork precedences (POST_OR)
            precMarker = 0;
            for ai = 1:sn.nacts
                probs = [];
                aidx = sn.ashift + ai;
                tidx = sn.parent(aidx);
                prob_ctr = 0;
                % for all successors
                for bidx=find(sn.graph(aidx,:))
                    if bidx > sn.ashift % ignore precedence between entries and activities
                        % Or pattern (POST_OR)
                        if full(sn.actposttype(bidx)) == ActivityPrecedenceType.POST_OR
                            if precMarker == 0 % start a new orjoin
                                precActs = java.util.ArrayList();
                                precMarker = aidx-sn.ashift;
                                precActs.add(sprintf("%s", sn.names{bidx}));
                                probs=full(sn.graph(aidx,bidx));
                            else
                                precActs.add(sprintf("%s", sn.names{bidx}));
                                probs(end+1)=full(sn.graph(aidx,bidx));
                            end
                        end
                    end
                end


                if precMarker > 0
                    probsMatrix = jline.util.matrix.Matrix(1,length(probs));
                    for i=1:length(probs)
                        probsMatrix.set(0,i-1,probs(i));
                    end
                    T{tidx-sn.tshift}.addPrecedence(jline.lang.layered.ActivityPrecedence.OrFork(sn.names{precMarker+sn.ashift}, precActs, probsMatrix));
                    precMarker = 0;
                end
            end

            %% AndFork precedences (POST_AND)
            precMarker = 0;
            postActs = '';
            for ai = 1:sn.nacts
                aidx = sn.ashift + ai;
                tidx = sn.parent(aidx);
                % for all successors
                for bidx=find(sn.graph(aidx,:))
                    if bidx > sn.ashift % ignore precedence between entries and activities
                        % Or pattern (POST_AND)
                        if full(sn.actposttype(bidx)) == ActivityPrecedenceType.POST_AND
                            if isempty(postActs)
                                postActs = java.util.ArrayList();
                                postActs.add(sprintf("%s", sn.names{bidx}));
                            else
                                postActs.add(sprintf("%s", sn.names{bidx}));
                            end

                            if precMarker == 0 % start a new orjoin
                                precMarker = aidx-sn.ashift;
                            end
                        end
                    end
                end
                if precMarker > 0
                    T{tidx-sn.tshift}.addPrecedence(jline.lang.layered.ActivityPrecedence.AndFork(sn.names{precMarker+sn.ashift}, postActs));
                    precMarker = 0;
                end
            end

            %% CacheAccess precedences (POST_CACHE)
            precMarker = 0;
            postActs = '';
            for ai = 1:sn.nacts
                aidx = sn.ashift + ai;
                tidx = sn.parent(aidx);
                % for all successors
                for bidx=find(sn.graph(aidx,:))
                    if bidx > sn.ashift % ignore precedence between entries and activities
                        % CacheAccess pattern (POST_CACHE)
                        if full(sn.actposttype(bidx)) == ActivityPrecedenceType.POST_CACHE
                            if isempty(postActs)
                                postActs = java.util.ArrayList();
                                postActs.add(sprintf("%s", sn.names{bidx}));
                            else
                                postActs.add(sprintf("%s", sn.names{bidx}));
                            end

                            if precMarker == 0 % start a new cache access
                                precMarker = aidx-sn.ashift;
                            end
                        end
                    end
                end
                if precMarker > 0
                    T{tidx-sn.tshift}.addPrecedence(jline.lang.layered.ActivityPrecedence.CacheAccess(sn.names{precMarker+sn.ashift}, postActs));
                    precMarker = 0;
                    postActs = '';
                end
            end

            %% OrJoin precedences (PRE_OR)
            precMarker = 0;
            for bi = sn.nacts:-1:1
                bidx = sn.ashift + bi;
                tidx = sn.parent(bidx);
                % for all predecessors
                for aidx=find(sn.graph(:,bidx))'
                    if aidx > sn.ashift % ignore precedence between entries and activities
                        % OrJoin pattern (PRE_OR)
                        if full(sn.actpretype(aidx)) == ActivityPrecedenceType.PRE_OR
                            if precMarker == 0 % start a new orjoin
                                precActs = java.util.ArrayList();
                                precMarker = bidx-sn.ashift;
                                precActs.add(sprintf("%s", sn.names{aidx}));
                            else
                                precActs.add(sprintf("%s", sn.names{aidx}));
                            end
                        end
                    end
                end
                if precMarker > 0
                    T{tidx-sn.tshift}.addPrecedence(jline.lang.layered.ActivityPrecedence.OrJoin(precActs, sn.names{precMarker+sn.ashift}));
                    precMarker = 0;
                end
            end

            %% AndJoin precedences (PRE_AND)
            precMarker = 0;
            for bi = sn.nacts:-1:1
                bidx = sn.ashift + bi;
                tidx = sn.parent(bidx);
                % for all predecessors
                for aidx=find(sn.graph(:,bidx))'
                    if aidx > sn.ashift % ignore precedence between entries and activities
                        % OrJoin pattern (PRE_AND)
                        if full(sn.actpretype(aidx)) == ActivityPrecedenceType.PRE_AND
                            if precMarker == 0 % start a new orjoin
                                precActs = java.util.ArrayList();
                                precMarker = bidx-sn.ashift;
                                precActs.add(sprintf("%s", sn.names{aidx}));
                            else
                                precActs.add(sprintf("%s", sn.names{aidx}));
                            end
                        end
                    end
                end
                if precMarker > 0
                    % Find quorum parameter from original precedence structure
                    postActName = sn.names{precMarker+sn.ashift};
                    localTaskIdx = tidx - sn.tshift;
                    quorum = [];
                    for ap = 1:length(line_layered_network.tasks{localTaskIdx}.precedences)
                        precedence = line_layered_network.tasks{localTaskIdx}.precedences(ap);
                        if precedence.preType == ActivityPrecedenceType.PRE_AND
                            % Check if this precedence contains our post activity
                            if any(strcmp(precedence.postActs, postActName))
                                quorum = precedence.preParams;
                                break;
                            end
                        end
                    end
                    
                    if isempty(quorum)
                        T{tidx-sn.tshift}.addPrecedence(jline.lang.layered.ActivityPrecedence.AndJoin(precActs, sn.names{precMarker+sn.ashift}));
                    else
                        T{tidx-sn.tshift}.addPrecedence(jline.lang.layered.ActivityPrecedence.AndJoin(precActs, sn.names{precMarker+sn.ashift}, quorum));
                    end
                    precMarker = 0;
                end
            end

        end

        function jdist = from_line_distribution(line_dist)
            if isa(line_dist, 'Exp')
                jdist = jline.lang.processes.Exp(line_dist.getParam(1).paramValue);
            elseif isa(line_dist, "APH")
                alpha = line_dist.getParam(1).paramValue;
                T = line_dist.getParam(2).paramValue;
                jline_alpha = java.util.ArrayList();
                for i = 1:length(alpha)
                    jline_alpha.add(alpha(i));
                end
                jline_T = JLINE.from_line_matrix(T);
                jdist = jline.lang.processes.APH(jline_alpha, jline_T);
            elseif isa(line_dist, 'Coxian')
                jline_mu = java.util.ArrayList();
                jline_phi = java.util.ArrayList();
                if length(line_dist.params) == 3
                    jline_mu.add(line_dist.getParam(1).paramValue);
                    jline_mu.add(line_dist.getParam(2).paramValue);
                    jline_phi.add(line_dist.getParam(3).paramValue);
                else
                    mu = line_dist.getParam(1).paramValue;
                    phi = line_dist.getParam(2).paramValue;
                    for i = 1:length(mu)
                        jline_mu.add(mu(i));
                    end
                    for i = 1:length(phi)
                        jline_phi.add(phi(i));
                    end
                end
                jdist = jline.lang.processes.Coxian(jline_mu, jline_phi);
            elseif isa(line_dist, 'Det')
                jdist = jline.lang.processes.Det(line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'DiscreteSampler')
                popularity_p = JLINE.from_line_matrix(line_dist.getParam(1).paramValue);
                popularity_val = JLINE.from_line_matrix(line_dist.getParam(2).paramValue);
                jdist = jline.lang.processes.DiscreteSampler(popularity_p, popularity_val);
            elseif isa(line_dist, 'Erlang')
                jdist = jline.lang.processes.Erlang(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Gamma')
                jdist = jline.lang.processes.Gamma(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, "HyperExp")
                jdist = jline.lang.processes.HyperExp(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue, line_dist.getParam(3).paramValue);
            elseif isa(line_dist, 'Lognormal')
                jdist = jline.lang.processes.Lognormal(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Pareto')
                jdist = jline.lang.processes.Pareto(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'MAP')
                D0 = line_dist.D(0);
                D1 = line_dist.D(1);
                jdist = jline.lang.processes.MAP(JLINE.from_line_matrix(D0), JLINE.from_line_matrix(D1));
            elseif isa(line_dist, 'MMPP2')
                lambda0 =  line_dist.getParam(1).paramValue;
                lambda1 =  line_dist.getParam(2).paramValue;
                sigma0 =  line_dist.getParam(3).paramValue;
                sigma1 =  line_dist.getParam(4).paramValue;
                jdist = jline.lang.processes.MMPP2(lambda0, lambda1, sigma0, sigma1);
            elseif isa(line_dist, 'PH')
                alpha = line_dist.getParam(1).paramValue;
                T = line_dist.getParam(2).paramValue;
                jdist = jline.lang.processes.PH(JLINE.from_line_matrix(alpha), JLINE.from_line_matrix(T));
            elseif isa(line_dist, 'Uniform')
                jdist = jline.lang.processes.Uniform(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Weibull')
                jdist = jline.lang.processes.Weibull(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Zipf')
                jdist = jline.lang.processes.Zipf(line_dist.getParam(3).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Immediate')
                jdist = jline.lang.processes.Immediate();
            elseif isempty(line_dist) || isa(line_dist, 'Disabled')
                jdist = jline.lang.processes.Disabled();
                return;
            elseif isa(line_dist, 'Replayer')
                jdist = jline.lang.processes.Replayer(line_dist.params{1}.paramValue);
            elseif isa(line_dist, 'Trace')
                jdist = jline.lang.processes.Trace(line_dist.params{1}.paramValue);
            else
                line_error(mfilename,'Distribution not supported by JLINE.');
            end
        end

        function matlab_dist = from_jline_distribution(jdist)
            if isa(jdist, 'jline.lang.processes.Exp')
                matlab_dist = Exp(jdist.getRate());
            elseif isa(jdist, 'jline.lang.processes.Det')
                matlab_dist = Det(jdist.getParam(1).getValue);
            elseif isa(jdist, 'jline.lang.processes.Erlang')
                matlab_dist = Erlang(jdist.getRate(),jdist.getNumberOfPhases());
            elseif isa(jdist, 'jline.lang.processes.Gamma')
                matlab_dist = Gamma(jdist.getParam(1).getValue,jdist.getParam(2).getValue);
            elseif isa(jdist, 'jline.lang.processes.HyperExp')
                matlab_dist = HyperExp(jdist.getParam(1).getValue, jdist.getParam(2).getValue, jdist.getParam(3).getValue);
            elseif isa(jdist, 'jline.lang.processes.Lognormal')
                matlab_dist = Lognormal(jdist.getParam(1).getValue,jdist.getParam(2).getValue);
            elseif isa(jdist, 'jline.lang.processes.Pareto')
                matlab_dist = Pareto(jdist.getParam(1).getValue,jdist.getParam(2).getValue);
            elseif isa(jdist, 'jline.lang.processes.Uniform')
                matlab_dist = Uniform(jdist.getParam(1).getValue,jdist.getParam(2).getValue);
            elseif isa(jdist, 'jline.lang.processes.Weibull')
                matlab_dist = Weibull(jdist.getParam(1).getValue,jdist.getParam(2).getValue);
            elseif isa(jdist, 'jline.lang.processes.APH')
                jalpha = jdist.getParam(1).getValue;
                alpha = zeros(1, jalpha.length);
                for i = 1:jalpha.length
                    alpha(i) = jalpha.get(i-1);
                end
                matlab_dist = APH(alpha, JLINE.from_jline_matrix(jdist.getParam(2).getValue));
            elseif isa(jdist, 'jline.lang.processes.PH')
                jalpha = jdist.getParam(1).getValue;
                alpha = zeros(1, jalpha.length);
                for i = 1:jalpha.length
                    alpha(i) = jalpha.get(i-1);
                end
                matlab_dist = PH(alpha, JLINE.from_jline_matrix(jdist.getParam(2).getValue));
            elseif isa(jdist, 'jline.lang.processes.Coxian')
                if jdist.getNumberOfPhases == 2
                    matlab_dist = Coxian([jdist.getParam(1).getValue.get(0), jdist.getParam(1).getValue.get(1)], [jdist.getParam(2).getValue.get(0),1]);
                else
                    jmu = jdist.getParam(1).getValue;
                    jphi = jdist.getParam(2).getValue;
                    mu = zeros(1, jmu.size);
                    phi = zeros(1, jphi.size);
                    for i = 1:jmu.size
                        mu(i) = jmu.get(i-1);
                    end
                    for i = 1:jphi.size
                        phi(i) = jphi.get(i-1);
                    end
                    matlab_dist = Coxian(mu, phi);
                end
            elseif isa(jdist, 'jline.lang.processes.Immediate')
                matlab_dist = Immediate();
            elseif isa(jdist, 'jline.lang.processes.Disabled')
                matlab_dist = Disabled();
            else
                line_error(mfilename,'Distribution not supported by JLINE.');
            end
        end

        function set_csMatrix(line_node, jnode, jclasses)
            nClasses = length(line_node.model.classes);
            csMatrix = jnode.initClassSwitchMatrix();
            for i = 1:nClasses
                for j = 1:nClasses
                    csMatrix.set(jclasses{i}, jclasses{j}, line_node.server.csFun(i,j,0,0));
                end
            end
            jnode.setClassSwitchingMatrix(csMatrix);
        end

        function set_service(line_node, jnode, job_classes)
            if (isa(line_node, 'Sink') || isa(line_node, 'Router') || isa(line_node, 'Cache') || isa(line_node, 'Logger') || isa(line_node, 'ClassSwitch') || isa(line_node, 'Fork') || isa(line_node, 'Join') || isa(line_node, 'Place') || isa(line_node, 'Transition'))
                return;
            end

            for n = 1 : length(job_classes)
                if (isa(line_node, 'Queue') || isa(line_node, 'Delay'))
                    matlab_dist = line_node.getService(job_classes{n});
                elseif (isa(line_node, 'Source'))
                    matlab_dist = line_node.getArrivalProcess(job_classes{n});
                else
                    line_error(mfilename,'Node not supported by JLINE.');
                end
                service_dist = JLINE.from_line_distribution(matlab_dist);

                if (isa(line_node,'Queue') || isa(line_node, 'Delay'))
                    jnode.setService(jnode.getModel().getClasses().get(n-1), service_dist, line_node.schedStrategyPar(n));
                elseif (isa(line_node, 'Source'))
                    jnode.setArrival(jnode.getModel().getClasses().get(n-1), service_dist);
                end
            end
        end

        function set_delayoff(line_node, jnode, job_classes)
            % Transfer setup and delayoff times from MATLAB Queue to Java Queue
            if ~isa(line_node, 'Queue')
                return;
            end

            % Check if setupTime property exists and is not empty
            if ~isprop(line_node, 'setupTime') || isempty(line_node.setupTime)
                return;
            end

            for n = 1 : length(job_classes)
                c = job_classes{n}.index;
                % Check if both setupTime and delayoffTime are set for this class
                if c <= length(line_node.setupTime) && ~isempty(line_node.setupTime{1, c}) && ...
                   c <= length(line_node.delayoffTime) && ~isempty(line_node.delayoffTime{1, c})
                    % Convert MATLAB distributions to Java distributions
                    setup_dist = JLINE.from_line_distribution(line_node.setupTime{1, c});
                    delayoff_dist = JLINE.from_line_distribution(line_node.delayoffTime{1, c});
                    % Set delayoff on the Java Queue
                    jnode.setDelayOff(jnode.getModel().getClasses().get(n-1), setup_dist, delayoff_dist);
                end
            end
        end

        function set_line_service(jline_node, line_node, job_classes, line_classes)
            if (isa(line_node,'Sink')) || isa(line_node, 'ClassSwitch')
                return;
            end
            for n = 1:job_classes.size()
                if (isa(line_node, 'Queue') || isa(line_node, 'Delay'))
                    jdist = jline_node.getServiceProcess(job_classes.get(n-1));
                    matlab_dist = JLINE.from_jline_distribution(jdist);
                    line_node.setService(line_classes{n}, matlab_dist);
                elseif (isa(line_node, 'Source'))
                    jdist = jline_node.getArrivalProcess(job_classes.get(n-1));
                    matlab_dist = JLINE.from_jline_distribution(jdist);
                    line_node.setArrival(line_classes{n}, matlab_dist);
                elseif (isa(line_node, 'Router'))
                    % no-op
                else
                    line_error(mfilename,'Node not supported by JLINE.');
                end
            end
        end

        function node_object = from_line_node(line_node, jnetwork, ~, forkNode, sn)
            % Handle optional sn argument
            if nargin < 5
                sn = [];
            end
            if isa(line_node, 'Delay')
                node_object = jline.lang.nodes.Delay(jnetwork, line_node.getName);
            elseif isa(line_node, 'Queue')
                switch line_node.schedStrategy
                    case SchedStrategy.INF
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.INF);
                    case SchedStrategy.FCFS
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFS);
                    case SchedStrategy.LCFS
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFS);
                    case SchedStrategy.SIRO
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.SIRO);
                    case SchedStrategy.SJF
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.SJF);
                    case SchedStrategy.LJF
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LJF);
                    case SchedStrategy.PS
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.PS);
                    case SchedStrategy.DPS
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.DPS);
                    case SchedStrategy.GPS
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.GPS);
                    case SchedStrategy.SEPT
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.SEPT);
                    case SchedStrategy.LEPT
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LEPT);
                    case {SchedStrategy.HOL, SchedStrategy.FCFSPRIO}
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFSPRIO);
                    case SchedStrategy.FORK
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FORK);
                    case SchedStrategy.EXT
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.EXT);
                    case SchedStrategy.REF
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.REF);
                    case SchedStrategy.LCFSPR
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFSPR);
                    case SchedStrategy.LCFSPI
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFSPI);
                    case SchedStrategy.LCFSPRIO
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFSPRIO);
                    case SchedStrategy.LCFSPRPRIO
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFSPRPRIO);
                    case SchedStrategy.LCFSPIPRIO
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFSPIPRIO);
                    case SchedStrategy.FCFSPR
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFSPR);
                    case SchedStrategy.FCFSPI
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFSPI);
                    case SchedStrategy.FCFSPRPRIO
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFSPRPRIO);
                    case SchedStrategy.FCFSPIPRIO
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFSPIPRIO);
                    case SchedStrategy.PSPRIO
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.PSPRIO);
                    case SchedStrategy.DPSPRIO
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.DPSPRIO);
                    case SchedStrategy.GPSPRIO
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.GPSPRIO);
                    case SchedStrategy.POLLING
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.POLLING);
                    case SchedStrategy.SRPT
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.SRPT);
                    case SchedStrategy.SRPTPRIO
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.SRPTPRIO);
                    case SchedStrategy.PSJF
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.PSJF);
                    case SchedStrategy.FB
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FB);
                    case SchedStrategy.LRPT
                        node_object = jline.lang.nodes.Queue(jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LRPT);
                end
                nservers = line_node.getNumberOfServers;
                if isinf(nservers)
                    node_object.setNumberOfServers(java.lang.Integer.MAX_VALUE);
                elseif nservers > 1
                    node_object.setNumberOfServers(line_node.getNumberOfServers);
                end
                if ~isempty(line_node.lldScaling)
                    node_object.setLoadDependence(JLINE.from_line_matrix(line_node.lldScaling));
                end
                if ~isempty(line_node.lcdScaling)
                    if isempty(sn)
                        line_error(mfilename, "Class-dependent models require sn struct for MATLAB-to-JAVA translation.");
                    end
                    node_object.setLimitedClassDependence(JLINE.handle_to_serializablefun(line_node.lcdScaling, sn));
                end
                % Transfer LJD (Limited Joint Dependence) scaling
                if ~isempty(line_node.ljdScaling) && ~isempty(line_node.ljdCutoffs)
                    jScaling = JLINE.from_line_matrix(line_node.ljdScaling(:));
                    jCutoffs = JLINE.from_line_matrix(line_node.ljdCutoffs(:));
                    node_object.setLimitedJointDependence(jScaling, jCutoffs);
                end
                % NOTE: LJCD (Limited Joint Class Dependence) is transferred in
                % from_line_network after classes are created
                % Set queue capacity if finite
                if ~isinf(line_node.cap)
                    node_object.setCapacity(line_node.cap);
                end
            elseif isa(line_node, 'Source')
                node_object = jline.lang.nodes.Source(jnetwork, line_node.getName);
            elseif isa(line_node, 'Sink')
                node_object = jline.lang.nodes.Sink(jnetwork, line_node.getName);
            elseif isa(line_node, 'Router')
                node_object = jline.lang.nodes.Router(jnetwork, line_node.getName);
            elseif isa(line_node, 'ClassSwitch')
                node_object = jline.lang.nodes.ClassSwitch(jnetwork, line_node.getName);
            elseif isa(line_node, 'Fork')
                node_object = jline.lang.nodes.Fork(jnetwork, line_node.name);
                node_object.setTasksPerLink(line_node.output.tasksPerLink);
            elseif isa(line_node, 'Join')
                node_object = jline.lang.nodes.Join(jnetwork, line_node.name, forkNode);
            elseif isa(line_node, 'Logger')
                node_object = jline.lang.nodes.Logger(jnetwork, line_node.name, [line_node.filePath,line_node.fileName]);
            elseif isa(line_node, 'Cache')
                nitems = line_node.items.nitems;
                switch line_node.replacestrategy
                    case ReplacementStrategy.RR
                        repStrategy = jline.lang.constant.ReplacementStrategy.RR;
                    case ReplacementStrategy.FIFO
                        repStrategy = jline.lang.constant.ReplacementStrategy.FIFO;
                    case ReplacementStrategy.SFIFO
                        repStrategy = jline.lang.constant.ReplacementStrategy.SFIFO;
                    case ReplacementStrategy.LRU
                        repStrategy = jline.lang.constant.ReplacementStrategy.LRU;
                end
                if ~isempty(line_node.graph)
                    node_object = jline.lang.nodes.Cache(jnetwork, line_node.name, nitems, JLINE.from_line_matrix(line_node.itemLevelCap), repStrategy, JLINE.from_line_matrix(line_node.graph));
                else
                    node_object = jline.lang.nodes.Cache(jnetwork, line_node.name, nitems, JLINE.from_line_matrix(line_node.itemLevelCap), repStrategy);
                end
            elseif isa(line_node, 'Place')
                node_object = jline.lang.nodes.Place(jnetwork, line_node.getName);
            elseif isa(line_node, 'Transition')
                node_object = jline.lang.nodes.Transition(jnetwork, line_node.getName);
                % Modes are added later in from_line_network after classes are created
            else
                line_error(mfilename,'Node not supported by JLINE.');
            end
        end

        function node_object = from_jline_node(jline_node, model, job_classes)
            if isa(jline_node, 'jline.lang.nodes.Delay')
                node_object = Delay(model, jline_node.getName.toCharArray');
            elseif isa(jline_node, 'jline.lang.nodes.Queue')
                schedStrategy = jline_node.getSchedStrategy;
                switch schedStrategy.name().toCharArray'
                    case 'INF'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.INF);
                    case 'FCFS'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.FCFS);
                    case 'LCFS'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.LCFS);
                    case 'SIRO'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.SIRO);
                    case 'SJF'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.SJF);
                    case 'LJF'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.LJF);
                    case 'PS'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.PS);
                    case 'DPS'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.DPS);
                    case 'GPS'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.GPS);
                    case 'SEPT'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.SEPT);
                    case 'LEPT'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.LEPT);
                    case {'HOL', 'FCFSPRIO'}
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.FCFSPRIO);
                    case 'FORK'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.FORK);
                    case 'EXT'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.EXT);
                    case 'REF'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.REF);
                    case 'LCFSPR'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.LCFSPR);
                    case 'SRPT'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.SRPT);
                    case 'SRPTPRIO'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.SRPTPRIO);
                    case 'PSJF'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.PSJF);
                    case 'FB'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.FB);
                    case 'LRPT'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.LRPT);
                end
                node_object.setNumberOfServers(jline_node.getNumberOfServers);
                if ~isempty(JLINE.from_jline_matrix(jline_node.getLimitedLoadDependence))
                    node_object.setLoadDependence(JLINE.from_jline_matrix(jline_node.getLimitedLoadDependence));
                end
            elseif isa(jline_node, 'jline.lang.nodes.Source')
                node_object = Source(model, jline_node.getName.toCharArray');
            elseif isa(jline_node, 'jline.lang.nodes.Sink')
                node_object = Sink(model, jline_node.getName.toCharArray');
            elseif isa(jline_node, 'jline.lang.nodes.Router')
                node_object = Router(model, jline_node.getName.toCharArray');
            elseif isa(jline_node, 'jline.lang.nodes.ClassSwitch')
                nClasses = job_classes.size;
                csMatrix = zeros(nClasses, nClasses);
                for r = 1:nClasses
                    for s = 1:nClasses
                        csMatrix(r,s) = jline_node.getServer.applyCsFun(r-1,s-1);
                    end
                end
                node_object = ClassSwitch(model, jline_node.getName.toCharArray', csMatrix);
            elseif isa(jline_node, 'jline.lang.nodes.Place')
                node_object = Place(model, jline_node.getName.toCharArray');
            elseif isa(jline_node, 'jline.lang.nodes.Transition')
                node_object = Transition(model, jline_node.getName.toCharArray');
                % Note: Mode configurations need to be set after classes are created
            else
                line_error(mfilename,'Node not supported by JLINE.');
            end
        end

        function node_class = from_line_class(line_class, jnetwork)
            % Check signal classes first (before their base classes)
            if isa(line_class, 'ClosedSignal')
                % ClosedSignal -> jline.lang.ClosedSignal
                jSignalType = jline.lang.constant.SignalType.fromID(line_class.signalType);
                node_class = jline.lang.ClosedSignal(jnetwork, line_class.getName, jSignalType, jnetwork.getNodeByName(line_class.refstat.getName), line_class.priority);
            elseif isa(line_class, 'Signal') || isa(line_class, 'OpenSignal')
                % Signal/OpenSignal -> jline.lang.Signal (includes CATASTROPHE type)
                jSignalType = jline.lang.constant.SignalType.fromID(line_class.signalType);
                node_class = jline.lang.Signal(jnetwork, line_class.getName, jSignalType, line_class.priority);
            elseif isa(line_class, 'OpenClass')
                node_class = jline.lang.OpenClass(jnetwork, line_class.getName, line_class.priority);
            elseif isa(line_class, 'SelfLoopingClass')
                node_class = jline.lang.SelfLoopingClass(jnetwork, line_class.getName, line_class.population, jnetwork.getNodeByName(line_class.refstat.getName), line_class.priority);
            elseif isa(line_class, 'ClosedClass')
                node_class = jline.lang.ClosedClass(jnetwork, line_class.getName, line_class.population, jnetwork.getNodeByName(line_class.refstat.getName), line_class.priority);
            else
                line_error(mfilename,'Class type not supported by JLINE.');
            end
            if line_class.isReferenceClass()
                node_class.setReferenceClass(true);
            end
        end

        function node_class = from_jline_class(jclass, model)
            if isa(jclass, 'jline.lang.OpenClass')
                node_class = OpenClass(model, jclass.getName.toCharArray', jclass.getPriority);
            elseif isa(jclass, 'jline.lang.SelfLoopingClass')
                node_class = SelfLoopingClass(model, jclass.getName.toCharArray', jclass.getNumberOfJobs, model.getNodeByName(jclass.getReferenceStation.getName), jclass.getPriority);
            elseif isa(jclass, 'jline.lang.ClosedClass')
                node_class = ClosedClass(model, jclass.getName.toCharArray', jclass.getNumberOfJobs, model.getNodeByName(jclass.getReferenceStation.getName), jclass.getPriority);
            else
                line_error(mfilename,'Class type not supported by JLINE.');
            end
        end

        function from_line_links(model, jmodel)
            connections = model.getConnectionMatrix();
            [m, ~] = size(connections);
            jnodes = jmodel.getNodes();
            jclasses = jmodel.getClasses();
            njclasses = jclasses.size();
            line_nodes = model.getNodes;
            sn = model.getStruct;

            % Build mapping from MATLAB node index to Java node index
            % (accounting for skipped auto-added ClassSwitch nodes)
            matlab2java_node_idx = zeros(1, length(line_nodes));
            jidx = 0;
            for i = 1:length(line_nodes)
                if isa(line_nodes{i}, 'ClassSwitch') && line_nodes{i}.autoAdded
                    matlab2java_node_idx(i) = -1;  % Mark as skipped
                else
                    matlab2java_node_idx(i) = jidx;
                    jidx = jidx + 1;
                end
            end

            %nodevisits = cellsum(sn.nodevisits);
            % [ ] Update to consider different weights/routing for classes
            if isempty(sn.rtorig)
                useLinkMethod = false; % this model did not call link()
            else
                jrt_matrix = jmodel.initRoutingMatrix();
                useLinkMethod = true;
            end

            % For models with auto-added ClassSwitch nodes, use sn.rtorig directly
            % to set up routing with proper class switching
            hasAutoCS = false;
            for i = 1:length(line_nodes)
                if isa(line_nodes{i}, 'ClassSwitch') && line_nodes{i}.autoAdded
                    hasAutoCS = true;
                    break;
                end
            end

            if useLinkMethod && hasAutoCS
                % Use sn.rtorig directly - it contains the full routing with class switching
                % sn.rtorig already excludes auto-added ClassSwitch nodes (it's based on nstations)
                % So we iterate over the rtorig matrix dimensions directly
                for r = 1:njclasses
                    for s = 1:njclasses
                        if ~isempty(sn.rtorig{r,s})
                            Prs = sn.rtorig{r,s};
                            [nrows, ncols] = size(Prs);
                            for i = 1:nrows
                                for j = 1:ncols
                                    if Prs(i,j) > 0
                                        % sn.rtorig uses station indices which map to non-CS nodes
                                        % Find the java node indices by matching station index to node
                                        jsrc_idx = i - 1;  % Direct mapping since rtorig excludes CS
                                        jdest_idx = j - 1;
                                        jrt_matrix.set(jclasses.get(r-1), jclasses.get(s-1), jnodes.get(jsrc_idx), jnodes.get(jdest_idx), Prs(i,j));
                                    end
                                end
                            end
                        end
                    end
                end
            else
                % Original logic for models without auto-added ClassSwitch
                for i = 1:m
                    line_node = line_nodes{i};

                    % Skip auto-added ClassSwitch nodes - Java will add them automatically
                    if isa(line_node, 'ClassSwitch') && line_node.autoAdded
                        continue;
                    end

                    jnode_idx = matlab2java_node_idx(i);
                    for k = 1:njclasses
                        output_strat = line_node.output.outputStrategy{k};
                        switch RoutingStrategy.fromText(output_strat{2})
                            case RoutingStrategy.DISABLED
                                jnodes.get(jnode_idx).setRouting(jclasses.get(k-1),jline.lang.constant.RoutingStrategy.DISABLED);
                            case RoutingStrategy.RAND
                                jnodes.get(jnode_idx).setRouting(jclasses.get(k-1),jline.lang.constant.RoutingStrategy.RAND);
                                outlinks_i=find(connections(i,:));
                                if ~useLinkMethod
                                    for j= outlinks_i(:)'
                                        jdest_idx = matlab2java_node_idx(j);
                                        if jdest_idx >= 0
                                            jmodel.addLink(jnodes.get(jnode_idx), jnodes.get(jdest_idx));
                                        end
                                    end
                                end
                            case RoutingStrategy.RROBIN
                                jnodes.get(jnode_idx).setRouting(jclasses.get(k-1),jline.lang.constant.RoutingStrategy.RROBIN);
                                outlinks_i=find(connections(i,:))';
                                if useLinkMethod
                                    line_error(mfilename,'RROBIN cannot be used together with the link() command.');
                                end
                                for j= outlinks_i(:)'
                                    jdest_idx = matlab2java_node_idx(j);
                                    if jdest_idx >= 0
                                        jmodel.addLink(jnodes.get(jnode_idx), jnodes.get(jdest_idx));
                                    end
                                end
                            case RoutingStrategy.WRROBIN
                                outlinks_i=find(connections(i,:))';
                                for j= outlinks_i(:)'
                                    jdest_idx = matlab2java_node_idx(j);
                                    if jdest_idx >= 0
                                        jmodel.addLink(jnodes.get(jnode_idx), jnodes.get(jdest_idx));
                                    end
                                end
                                if useLinkMethod
                                    line_error(mfilename,'RROBIN cannot be used together with the link() command.');
                                end
                                for j= 1:length(output_strat{3})
                                    node_target = jmodel.getNodeByName(output_strat{3}{j}{1}.getName());
                                    weight = output_strat{3}{j}{2};
                                    jnodes.get(jnode_idx).setRouting(jclasses.get(k-1),jline.lang.constant.RoutingStrategy.WRROBIN, node_target, weight);
                                end
                            case RoutingStrategy.PROB
                                outlinks_i=find(connections(i,:));
                                jnodes.get(jnode_idx).setRouting(jclasses.get(k-1), jline.lang.constant.RoutingStrategy.PROB);
                                if ~useLinkMethod
                                    for j= outlinks_i(:)'
                                        jdest_idx = matlab2java_node_idx(j);
                                        if jdest_idx >= 0
                                            jmodel.addLink(jnodes.get(jnode_idx), jnodes.get(jdest_idx));
                                        end
                                    end
                                end
                                if length(output_strat) >= 3
                                    probabilities = output_strat{3};
                                    for j = 1:length(probabilities)
                                        dest_idx = probabilities{j}{1}.index;
                                        jdest_idx = matlab2java_node_idx(dest_idx);
                                        if (connections(i, dest_idx) ~= 0) && jdest_idx >= 0
                                            if useLinkMethod
                                                jrt_matrix.set(jclasses.get(k-1), jclasses.get(k-1), jnodes.get(jnode_idx), jnodes.get(jdest_idx), probabilities{j}{2});
                                            else
                                                jnodes.get(jnode_idx).setProbRouting(jclasses.get(k-1), jnodes.get(jdest_idx), probabilities{j}{2});
                                            end
                                        end
                                    end
                                end
                            otherwise
                                line_warning(mfilename, sprintf('''%s'' routing strategy not supported by JLINE, setting as Disabled.\n',output_strat{2}));
                                jnodes.get(jnode_idx).setRouting(jclasses.get(k-1),jline.lang.constant.RoutingStrategy.DISABLED);
                        end
                    end
                end
            end
            if useLinkMethod
                jmodel.link(jrt_matrix);
                % Align the sn.rtorig be the same, treating artificial
                % ClassSwitch nodes as if they were explicitly specified
                jsn = jmodel.getStruct(true);
                rtorig = java.util.HashMap();
                if ~isempty(model.sn.rtorig)
                    if iscell(model.sn.rtorig)
                        for r = 1:njclasses
                            sub_rtorig = java.util.HashMap();
                            for s = 1:njclasses
                                sub_rtorig.put(jclasses.get(s-1), JLINE.from_line_matrix(model.sn.rtorig{r,s}));
                            end
                            rtorig.put(jclasses.get(r-1), sub_rtorig);
                        end
                    end
                end
                jsn.rtorig = rtorig;
            end
        end

        function model = from_jline_routing(model, jnetwork)
            jnodes = jnetwork.getNodes();
            jclasses = jnetwork.getClasses();
            %n_classes = jclasses.size();
            n_nodes = jnodes.size();
            network_nodes = model.getNodes;
            network_classes = model.getClasses;

            connections = JLINE.from_jline_matrix(jnetwork.getConnectionMatrix());
            [row,col] = find(connections);
            for i=1:length(row)
                model.addLink(network_nodes{row(i)},network_nodes{col(i)});
            end

            for n = 1 : n_nodes
                jnode = jnodes.get(n-1);
                output_strategies = jnode.getOutputStrategies();
                n_strategies = output_strategies.size();
                for m = 1 : n_strategies
                    output_strat = output_strategies.get(m-1);
                    routing_strat = output_strat.getRoutingStrategy;
                    routing_strat_classidx = output_strat.getJobClass.getIndex();
                    switch char(routing_strat)
                        case 'RAND'
                            network_nodes{n}.setRouting(network_classes{routing_strat_classidx}, RoutingStrategy.RAND);
                        case 'RROBIN'
                            network_nodes{n}.setRouting(network_classes{routing_strat_classidx}, RoutingStrategy.RROBIN);
                        case 'DISABLED'
                            network_nodes{n}.setRouting(network_classes{routing_strat_classidx}, RoutingStrategy.RROBIN);
                    end
                end
            end
        end

        function model = from_jline_links(model, jnetwork)
            P = model.initRoutingMatrix;
            jnodes = jnetwork.getNodes();
            jclasses = jnetwork.getClasses();
            n_classes = jclasses.size();
            n_nodes = jnodes.size();

            for n = 1 : n_nodes
                jnode = jnodes.get(n-1);
                output_strategies = jnode.getOutputStrategies();
                n_strategies = output_strategies.size();
                for m = 1 : n_strategies
                    output_strat = output_strategies.get(m-1);
                    dest = output_strat.getDestination();
                    if~isempty(dest) % disabled strategy
                        in_idx = jnetwork.getNodeIndex(jnode)+1;
                        out_idx = jnetwork.getNodeIndex(dest)+1;
                        if n_classes == 1
                            P{1}(in_idx,out_idx) = output_strat.getProbability();
                        else
                            strat_class = output_strat.getJobClass();
                            class_idx = jnetwork.getJobClassIndex(strat_class)+1;
                            P{class_idx,class_idx}(in_idx,out_idx) = output_strat.getProbability();
                        end
                    end
                end
            end
            model.link(P);

            %Align the sn.rtorig be the same (Assume Java network is
            %created by calling Network.link)
            sn = jnetwork.getStruct;
            rtorig = cell(n_classes, n_classes);
            for r = 1:n_classes
                for s = 1:n_classes
                    rtorig{r,s} = JLINE.from_jline_matrix(sn.rtorig.get(jclasses.get(r-1)).get(jclasses.get(s-1)));
                end
            end
            model.sn.rtorig = rtorig;
        end

        function [jnetwork] = from_line_network(model)
            %w = warning;
            %warning('off');
            sn = model.getStruct;

            jnetwork = jline.lang.Network(model.getName);
            line_nodes = model.getNodes;
            line_classes = model.getClasses;

            jnodes = cell(1,length(line_nodes));
            jclasses = cell(1,length(line_classes));

            for n = 1 : length(line_nodes)
                % Skip auto-added ClassSwitch nodes - Java's link() will add them automatically
                if isa(line_nodes{n}, 'ClassSwitch') && line_nodes{n}.autoAdded
                    continue;
                end
                if isa(line_nodes{n}, 'Join')
                    jnodes{n} = JLINE.from_line_node(line_nodes{n}, jnetwork, line_classes, jnodes{line_nodes{n}.joinOf.index}, sn);
                else
                    jnodes{n} = JLINE.from_line_node(line_nodes{n}, jnetwork, line_classes, [], sn);
                end
            end

            for n = 1 : length(line_classes)
                jclasses{n} = JLINE.from_line_class(line_classes{n}, jnetwork);
            end

            % Set up forJobClass associations for signal classes
            for n = 1 : length(line_classes)
                if (isa(line_classes{n}, 'Signal') || isa(line_classes{n}, 'OpenSignal') || isa(line_classes{n}, 'ClosedSignal'))
                    if ~isempty(line_classes{n}.targetJobClass)
                        targetIdx = line_classes{n}.targetJobClass.index;
                        jclasses{n}.forJobClass(jclasses{targetIdx});
                    end
                end
            end

            for n = 1: length(jnodes)
                if isempty(jnodes{n})
                    continue;  % Skip nodes that were not converted (e.g., auto-added ClassSwitch)
                end
                JLINE.set_service(line_nodes{n}, jnodes{n}, line_classes);
                JLINE.set_delayoff(line_nodes{n}, jnodes{n}, line_classes);
            end

            % Set drop rules for stations (after classes are created)
            for n = 1: length(jnodes)
                if isempty(jnodes{n})
                    continue;  % Skip nodes that were not converted
                end
                if isa(line_nodes{n}, 'Station') && ~isa(line_nodes{n}, 'Source')
                    for r = 1:length(line_classes)
                        if length(line_nodes{n}.dropRule) >= r && ~isempty(line_nodes{n}.dropRule(r))
                            dropRule = line_nodes{n}.dropRule(r);
                            switch dropRule
                                case DropStrategy.DROP
                                    jnodes{n}.setDropRule(jclasses{r}, jline.lang.constant.DropStrategy.Drop);
                                case DropStrategy.BAS
                                    jnodes{n}.setDropRule(jclasses{r}, jline.lang.constant.DropStrategy.BlockingAfterService);
                                case DropStrategy.WAITQ
                                    jnodes{n}.setDropRule(jclasses{r}, jline.lang.constant.DropStrategy.WaitingQueue);
                            end
                        end
                    end
                end
            end

            % Transfer LJCD (Limited Joint Class Dependence) scaling for stations
            % This must be done after classes are created since LJCD is indexed per-class
            for n = 1:length(jnodes)
                if isempty(jnodes{n})
                    continue;
                end
                if isa(line_nodes{n}, 'Station') && ~isempty(line_nodes{n}.ljcdScaling) && ~isempty(line_nodes{n}.ljcdCutoffs)
                    jCutoffs = JLINE.from_line_matrix(line_nodes{n}.ljcdCutoffs(:));
                    K = length(line_nodes{n}.ljcdScaling);
                    jScalingMap = java.util.HashMap();
                    for c = 1:K
                        if ~isempty(line_nodes{n}.ljcdScaling{c})
                            jScalingVec = JLINE.from_line_matrix(line_nodes{n}.ljcdScaling{c}(:));
                            jScalingMap.put(jclasses{c}, jScalingVec);
                        end
                    end
                    jnodes{n}.setLimitedJointClassDependence(jScalingMap, jCutoffs);
                end
            end

            % Set polling type and switchover times for polling queues
            for n = 1: length(jnodes)
                if isempty(jnodes{n})
                    continue;
                end
                if isa(line_nodes{n}, 'Queue') && line_nodes{n}.schedStrategy == SchedStrategy.POLLING
                    % Set polling type
                    if ~isempty(line_nodes{n}.pollingType) && ~isempty(line_nodes{n}.pollingType{1})
                        pollingType = line_nodes{n}.pollingType{1};
                        switch pollingType
                            case PollingType.GATED
                                jPollingType = jline.lang.constant.PollingType.GATED;
                            case PollingType.EXHAUSTIVE
                                jPollingType = jline.lang.constant.PollingType.EXHAUSTIVE;
                            case PollingType.KLIMITED
                                jPollingType = jline.lang.constant.PollingType.KLIMITED;
                        end
                        if pollingType == PollingType.KLIMITED && ~isempty(line_nodes{n}.pollingPar)
                            jnodes{n}.setPollingType(jPollingType, int32(line_nodes{n}.pollingPar));
                        else
                            jnodes{n}.setPollingType(jPollingType);
                        end
                    end
                    % Set switchover times
                    if ~isempty(line_nodes{n}.switchoverTime)
                        for r = 1:length(line_classes)
                            if length(line_nodes{n}.switchoverTime) >= r && ~isempty(line_nodes{n}.switchoverTime{r})
                                soTime = line_nodes{n}.switchoverTime{r};
                                if ~isa(soTime, 'Immediate')
                                    jnodes{n}.setSwitchover(jclasses{r}, JLINE.from_line_distribution(soTime));
                                end
                            end
                        end
                    end
                end
            end

            for n = 1: length(jnodes)
                if isempty(jnodes{n})
                    continue;  % Skip nodes that were not converted
                end
                if isa(line_nodes{n},"ClassSwitch") && ~line_nodes{n}.autoAdded
                    % Only set csMatrix for user-defined ClassSwitch nodes (not auto-added)
                    JLINE.set_csMatrix(line_nodes{n}, jnodes{n}, jclasses);
                elseif isa(line_nodes{n},"Join")
                    jnodes{n}.initJoinJobClasses();
                elseif isa(line_nodes{n},"Cache")
                    for r = 1 : sn.nclasses
                        if length(line_nodes{n}.server.hitClass) >= r && ~isempty(line_nodes{n}.server.hitClass(r))
                            if ~isa(line_nodes{n}.popularity{r},'Disabled')
                                jnodes{n}.setRead(jclasses{r}, JLINE.from_line_distribution(line_nodes{n}.popularity{r}));
                                jnodes{n}.setHitClass(jclasses{r}, jclasses{line_nodes{n}.server.hitClass(r)});
                                jnodes{n}.setMissClass(jclasses{r}, jclasses{line_nodes{n}.server.missClass(r)});
                            end
                        end
                    end
                    % Transfer accessProb from MATLAB to Java
                    if ~isempty(line_nodes{n}.accessProb)
                        accessProbMat = line_nodes{n}.accessProb;
                        [K1, K2] = size(accessProbMat);
                        jAccessProb = javaArray('jline.util.matrix.Matrix', K1, K2);
                        for k1 = 1:K1
                            for k2 = 1:K2
                                if ~isempty(accessProbMat{k1, k2})
                                    jAccessProb(k1, k2) = JLINE.from_line_matrix(accessProbMat{k1, k2});
                                end
                            end
                        end
                        jnodes{n}.setAccessProb(jAccessProb);
                    end
                elseif isa(line_nodes{n}, "Transition")
                    % First, add modes (must be done after classes are created)
                    for m = 1:line_nodes{n}.getNumberOfModes()
                        modeName = line_nodes{n}.modeNames{m};
                        jmode = jnodes{n}.addMode(modeName);
                        % Set timing strategy
                        switch line_nodes{n}.timingStrategies(m)
                            case TimingStrategy.TIMED
                                jnodes{n}.setTimingStrategy(jmode, jline.lang.constant.TimingStrategy.TIMED);
                            case TimingStrategy.IMMEDIATE
                                jnodes{n}.setTimingStrategy(jmode, jline.lang.constant.TimingStrategy.IMMEDIATE);
                        end
                        % Set distribution
                        jnodes{n}.setDistribution(jmode, JLINE.from_line_distribution(line_nodes{n}.distributions{m}));
                        % Set firing weights and priorities
                        jnodes{n}.setFiringWeights(jmode, line_nodes{n}.firingWeights(m));
                        jnodes{n}.setFiringPriorities(jmode, int32(line_nodes{n}.firingPriorities(m)));
                        % Set number of servers
                        jnodes{n}.setNumberOfServers(jmode, java.lang.Integer(line_nodes{n}.numberOfServers(m)));
                    end
                    % Now set enabling conditions, inhibiting conditions, and firing outcomes
                    jmodes = jnodes{n}.getModes();
                    for m = 1:line_nodes{n}.getNumberOfModes()
                        jmode = jmodes.get(m-1);
                        enabCond = line_nodes{n}.enablingConditions{m};
                        inhibCond = line_nodes{n}.inhibitingConditions{m};
                        firingOut = line_nodes{n}.firingOutcomes{m};
                        for r = 1:sn.nclasses
                            for i = 1:length(line_nodes)
                                % Set enabling conditions
                                if enabCond(i, r) > 0 && isa(line_nodes{i}, 'Place')
                                    jnodes{n}.setEnablingConditions(jmode, jclasses{r}, jnodes{i}, enabCond(i, r));
                                end
                                % Set inhibiting conditions
                                if inhibCond(i, r) < Inf && isa(line_nodes{i}, 'Place')
                                    jnodes{n}.setInhibitingConditions(jmode, jclasses{r}, jnodes{i}, inhibCond(i, r));
                                end
                                % Set firing outcomes
                                if firingOut(i, r) ~= 0
                                    jnodes{n}.setFiringOutcome(jmode, jclasses{r}, jnodes{i}, firingOut(i, r));
                                end
                            end
                        end
                    end
                end
            end

            % Assume JLINE and LINE network are both created via link
            JLINE.from_line_links(model, jnetwork);

            % Transfer finite capacity regions from MATLAB to Java
            if ~isempty(model.regions)
                for f = 1:length(model.regions)
                    fcr = model.regions{f};
                    % Convert MATLAB node list to Java list
                    javaNodeList = java.util.ArrayList();
                    for i = 1:length(fcr.nodes)
                        matlabNode = fcr.nodes{i};
                        % Find corresponding Java node by name
                        nodeName = matlabNode.getName();
                        for j = 1:length(line_nodes)
                            if strcmp(line_nodes{j}.getName(), nodeName) && ~isempty(jnodes{j})
                                javaNodeList.add(jnodes{j});
                                break;
                            end
                        end
                    end
                    % Create Java FCR
                    jfcr = jnetwork.addRegion(javaNodeList);
                    % Set global max jobs
                    if fcr.globalMaxJobs > 0 && ~isinf(fcr.globalMaxJobs)
                        jfcr.setGlobalMaxJobs(fcr.globalMaxJobs);
                    end
                    % Set per-class max jobs and drop rules
                    for r = 1:length(line_classes)
                        if length(fcr.classMaxJobs) >= r && fcr.classMaxJobs(r) > 0 && ~isinf(fcr.classMaxJobs(r))
                            jfcr.setClassMaxJobs(jclasses{r}, fcr.classMaxJobs(r));
                        end
                        if length(fcr.dropRule) >= r
                            % Convert MATLAB DropStrategy numeric to Java DropStrategy enum
                            jDropStrategy = jline.lang.constant.DropStrategy.fromID(fcr.dropRule(r));
                            jfcr.setDropRule(jclasses{r}, jDropStrategy);
                        end
                    end
                end
            end

            jnetwork.initDefault;
            for n = 1: length(line_nodes)
                if isempty(jnodes{n})
                    continue;  % Skip nodes that were not converted
                end
                if line_nodes{n}.isStateful
                    jnodes{n}.setState(JLINE.from_line_matrix(line_nodes{n}.getState));
                    jnodes{n}.setStateSpace(JLINE.from_line_matrix(line_nodes{n}.getStateSpace));
                    jnodes{n}.setStatePrior(JLINE.from_line_matrix(line_nodes{n}.getStatePrior));
                end
            end
            % Force struct refresh so sn.state reflects updated node states
            jnetwork.setHasStruct(false);

        end

        function jnetwork = line_to_jline(model)
            jnetwork = LINE2JLINE(model);
        end

        function model = jline_to_line(jnetwork)
            if isa(jnetwork,'JNetwork')
                jnetwork = jnetwork.obj;
            end
            %javaaddpath(jar_loc);
            model = Network(char(jnetwork.getName));
            network_nodes = jnetwork.getNodes;
            job_classes = jnetwork.getClasses;

            line_nodes = cell(network_nodes.size,1);
            line_classes = cell(job_classes.size,1);


            for n = 1 : network_nodes.size
                if ~isa(network_nodes.get(n-1), 'jline.lang.nodes.ClassSwitch')
                    line_nodes{n} = JLINE.from_jline_node(network_nodes.get(n-1), model, job_classes);
                end
            end

            for n = 1 : job_classes.size
                line_classes{n} = JLINE.from_jline_class(job_classes.get(n-1), model);
            end

            for n = 1 : network_nodes.size
                if isa(network_nodes.get(n-1), 'jline.lang.nodes.ClassSwitch')
                    line_nodes{n} = JLINE.from_jline_node(network_nodes.get(n-1), model, job_classes);
                end
            end

            for n = 1 : network_nodes.size
                JLINE.set_line_service(network_nodes.get(n-1), line_nodes{n}, job_classes, line_classes);
            end

            if ~isempty(jnetwork.getStruct.rtorig)
                % Use link() method
                model = JLINE.from_jline_links(model, jnetwork);
            else
                % Do not use link() method
                model = JLINE.from_jline_routing(model, jnetwork);
            end
        end

        function matrix = arraylist_to_matrix(jline_matrix)
            if isempty(jline_matrix)
                matrix = [];
            else
                matrix = zeros(jline_matrix.size(), 1);
                for row = 1:jline_matrix.size()
                    matrix(row, 1) = jline_matrix.get(row-1);
                end
            end
        end

        function matrix = from_jline_matrix(jline_matrix)
            if isempty(jline_matrix)
                matrix = [];
            else
                matrix = zeros(jline_matrix.getNumRows(), jline_matrix.getNumCols());
                for row = 1:jline_matrix.getNumRows()
                    for col = 1:jline_matrix.getNumCols()
                        val = jline_matrix.get(row-1, col-1);
                        if (val >= 33333333 && val <= 33333334)
                            matrix(row, col) = GlobalConstants.Immediate;
                        elseif (val >= -33333334 && val <= -33333333)
                            matrix(row, col) = -GlobalConstants.Immediate;
                        elseif (val >= 2147483647 - 1) % Integer.MAX_VALUE with -1 tolerance
                            matrix(row, col) = Inf;
                        elseif (val <= -2147483648 + 1) % Integer.MIN_VALUE with +1 tolerance
                            matrix(row, col) = -Inf;
                        else
                            matrix(row, col) = val;
                        end
                    end
                end
            end
        end

        function jline_matrix = from_line_matrix(matrix)
            [rows, cols] = size(matrix);
            jline_matrix = jline.util.matrix.Matrix(rows, cols);
            for row = 1:rows
                for col = 1:cols
                    if matrix(row,col) ~= 0
                        jline_matrix.set(row-1, col-1, matrix(row, col));
                    end
                end
            end
        end

        function lsn = from_jline_struct_layered(jlayerednetwork, jlsn)
            lsn = LayeredNetworkStruct();
            lsn.nidx= jlsn.nidx;
            lsn.nhosts= jlsn.nhosts;
            lsn.ntasks= jlsn.ntasks;
            lsn.nentries= jlsn.nentries;
            lsn.nacts= jlsn.nacts;
            lsn.ncalls= jlsn.ncalls;
            lsn.hshift= jlsn.hshift;
            lsn.tshift= jlsn.tshift;
            lsn.eshift= jlsn.eshift;
            lsn.ashift= jlsn.ashift;
            lsn.cshift= jlsn.cshift;
            for h=1:jlsn.nhosts
                lsn.tasksof{h,1} = JLINE.arraylist_to_matrix(jlsn.tasksof.get(uint32(h)))';
            end
            for t=1:jlsn.ntasks
                lsn.entriesof{lsn.tshift+t,1} = JLINE.arraylist_to_matrix(jlsn.entriesof.get(uint32(jlsn.tshift+t)))';
            end
            for t=1:(jlsn.ntasks+jlsn.nentries)
                lsn.actsof{lsn.tshift+t,1} = JLINE.arraylist_to_matrix(jlsn.actsof.get(uint32(jlsn.tshift+t)))';
            end
            for a=1:jlsn.nacts
                lsn.callsof{lsn.ashift+a,1} = JLINE.arraylist_to_matrix(jlsn.callsof.get(uint32(jlsn.ashift+a)))';
            end
            for i = 1:jlsn.sched.size
                lsn.sched(i,1) = SchedStrategy.(char(jlsn.sched.get(uint32(i))));
            end
            for i = 1:jlsn.names.size
                lsn.names{i,1} = jlsn.names.get(uint32(i));
                lsn.hashnames{i,1} = jlsn.hashnames.get(uint32(i));
            end
            lsn.mult = JLINE.from_jline_matrix(jlsn.mult);
            lsn.mult = lsn.mult(2:(lsn.eshift+1))'; % remove 0-padding
            lsn.maxmult = JLINE.from_jline_matrix(jlsn.maxmult);
            lsn.maxmult = lsn.maxmult(2:(lsn.eshift+1))'; % remove 0-padding

            lsn.repl = JLINE.from_jline_matrix(jlsn.repl)';
            lsn.repl = lsn.repl(2:end); % remove 0-padding
            lsn.type = JLINE.from_jline_matrix(jlsn.type)';
            lsn.type = lsn.type(2:end); % remove 0-padding
            lsn.parent = JLINE.from_jline_matrix(jlsn.parent);
            lsn.parent = lsn.parent(2:end); % remove 0-padding
            lsn.nitems = JLINE.from_jline_matrix(jlsn.nitems);
            % Ensure proper column vector format matching MATLAB's (nhosts+ntasks+nentries) x 1
            if isrow(lsn.nitems)
                lsn.nitems = lsn.nitems(2:end)'; % remove 0-padding and transpose
            else
                lsn.nitems = lsn.nitems(2:end); % remove 0-padding (already column)
            end
            % Ensure correct size
            expectedSize = lsn.nhosts + lsn.ntasks + lsn.nentries;
            if length(lsn.nitems) < expectedSize
                lsn.nitems(expectedSize,1) = 0;
            elseif length(lsn.nitems) > expectedSize
                lsn.nitems = lsn.nitems(1:expectedSize);
            end
            lsn.replacestrat = JLINE.from_jline_matrix(jlsn.replacestrat);
            lsn.replacestrat = lsn.replacestrat(2:end)'; % remove 0-padding
            for i = 1:jlsn.callnames.size
                lsn.callnames{i,1} = jlsn.callnames.get(uint32(i));
                lsn.callhashnames{i,1} = jlsn.callhashnames.get(uint32(i));
            end
            for i = 1:jlsn.calltype.size % calltype may be made into a matrix in Java
                ct = char(jlsn.calltype.get(uint32(i)));
                lsn.calltype(i) = CallType.(ct);
            end
            lsn.calltype = sparse(lsn.calltype'); % remove 0-padding
            lsn.callpair = JLINE.from_jline_matrix(jlsn.callpair);
            lsn.callpair = lsn.callpair(2:end,2:end); % remove 0-paddings
            if isempty(lsn.callpair)
                lsn.callpair=[];
            end
            lsn.actpretype = sparse(JLINE.from_jline_matrix(jlsn.actpretype)');
            lsn.actpretype = lsn.actpretype(2:end); % remove 0-padding
            lsn.actposttype = sparse(JLINE.from_jline_matrix(jlsn.actposttype)');
            lsn.actposttype = lsn.actposttype(2:end); % remove 0-padding
            lsn.graph = JLINE.from_jline_matrix(jlsn.graph);
            lsn.graph = lsn.graph(2:end,2:end); % remove 0-paddings
            lsn.dag = JLINE.from_jline_matrix(jlsn.dag);
            lsn.dag = lsn.dag(2:end,2:end); % remove 0-paddings
            lsn.taskgraph = JLINE.from_jline_matrix(jlsn.taskgraph);
            lsn.taskgraph = sparse(lsn.taskgraph(2:end,2:end)); % remove 0-paddings
            lsn.replygraph = JLINE.from_jline_matrix(jlsn.replygraph);
            lsn.replygraph = logical(lsn.replygraph(2:end,2:end)); % remove 0-paddings
            lsn.iscache = JLINE.from_jline_matrix(jlsn.iscache);
            % Ensure proper column vector format matching MATLAB's (nhosts+ntasks) x 1
            expectedCacheSize = lsn.nhosts + lsn.ntasks;
            if isrow(lsn.iscache)
                if length(lsn.iscache) > expectedCacheSize
                    lsn.iscache = lsn.iscache(2:(expectedCacheSize+1))'; % remove 0-padding and transpose
                else
                    lsn.iscache = lsn.iscache'; % just transpose
                end
            end
            % Ensure correct size
            if length(lsn.iscache) < expectedCacheSize
                lsn.iscache(expectedCacheSize,1) = 0;
            elseif length(lsn.iscache) > expectedCacheSize
                lsn.iscache = lsn.iscache(1:expectedCacheSize);
            end
            lsn.iscaller = JLINE.from_jline_matrix(jlsn.iscaller);
            lsn.iscaller = full(lsn.iscaller(2:end,2:end)); % remove 0-paddings
            lsn.issynccaller = JLINE.from_jline_matrix(jlsn.issynccaller);
            lsn.issynccaller = full(lsn.issynccaller(2:end,2:end)); % remove 0-paddings
            lsn.isasynccaller = JLINE.from_jline_matrix(jlsn.isasynccaller);
            lsn.isasynccaller = full(lsn.isasynccaller(2:end,2:end)); % remove 0-paddings
            lsn.isref = JLINE.from_jline_matrix(jlsn.isref);
            lsn.isref = lsn.isref(2:end)'; % remove 0-paddings
        end

        function sn = from_jline_struct(jnetwork, jsn)
            %lst and rtfun are not implemented
            %Due to the transformation of Java lambda to matlab function
            if nargin<2
                jsn = jnetwork.getStruct(false);
            end
            jclasses = jnetwork.getClasses();
            jnodes = jnetwork.getNodes();
            jstateful = jnetwork.getStatefulNodes();
            jstations = jnetwork.getStations();
            sn = NetworkStruct();

            sn.nnodes = jsn.nnodes;
            sn.nclasses = jsn.nclasses;
            sn.nclosedjobs = jsn.nclosedjobs;
            sn.nstations = jsn.nstations;
            sn.nstateful = jsn.nstateful;
            sn.nchains = jsn.nchains;

            sn.refstat = JLINE.from_jline_matrix(jsn.refstat) + 1;
            sn.njobs = JLINE.from_jline_matrix(jsn.njobs);
            sn.nservers = JLINE.from_jline_matrix(jsn.nservers);
            sn.connmatrix = JLINE.from_jline_matrix(jsn.connmatrix);
            % Fix for Java getConnectionMatrix bug: ensure connmatrix is nnodes x nnodes
            if size(sn.connmatrix,1) < sn.nnodes
                sn.connmatrix(sn.nnodes,1) = 0;
            end
            if size(sn.connmatrix,2) < sn.nnodes
                sn.connmatrix(1,sn.nnodes) = 0;
            end
            sn.scv = JLINE.from_jline_matrix(jsn.scv);
            sn.isstation = logical(JLINE.from_jline_matrix(jsn.isstation));
            sn.isstateful = logical(JLINE.from_jline_matrix(jsn.isstateful));
            sn.isstatedep = logical(JLINE.from_jline_matrix(jsn.isstatedep));
            sn.nodeToStateful = JLINE.from_jline_matrix(jsn.nodeToStateful)+1;
            sn.nodeToStateful(sn.nodeToStateful==0) = nan;
            sn.nodeToStation = JLINE.from_jline_matrix(jsn.nodeToStation)+1;
            sn.nodeToStation(sn.nodeToStation==0) = nan;
            sn.stationToNode = JLINE.from_jline_matrix(jsn.stationToNode)+1;
            sn.stationToNode(sn.stationToNode==0) = nan;
            sn.stationToStateful = JLINE.from_jline_matrix(jsn.stationToStateful)+1;
            sn.stationToStateful(sn.stationToStateful==0) = nan;
            sn.statefulToStation = JLINE.from_jline_matrix(jsn.statefulToStation)+1;
            sn.statefulToStation(sn.statefulToStation==0) = nan;
            sn.statefulToNode = JLINE.from_jline_matrix(jsn.statefulToNode)+1;
            sn.statefulToNode(sn.statefulToNode==0) = nan;
            sn.rates = JLINE.from_jline_matrix(jsn.rates);
            sn.fj = JLINE.from_jline_matrix(jsn.fj);
            sn.classprio = JLINE.from_jline_matrix(jsn.classprio);
            sn.phases = JLINE.from_jline_matrix(jsn.phases);
            sn.phasessz = JLINE.from_jline_matrix(jsn.phasessz);
            sn.phaseshift = JLINE.from_jline_matrix(jsn.phaseshift);
            sn.schedparam = JLINE.from_jline_matrix(jsn.schedparam);
            sn.chains = logical(JLINE.from_jline_matrix(jsn.chains));
            sn.rt = JLINE.from_jline_matrix(jsn.rt);
            sn.nvars = JLINE.from_jline_matrix(jsn.nvars);
            sn.rtnodes = JLINE.from_jline_matrix(jsn.rtnodes);
            sn.csmask = logical(JLINE.from_jline_matrix(jsn.csmask));
            sn.isslc = logical(JLINE.from_jline_matrix(jsn.isslc));
            sn.cap = JLINE.from_jline_matrix(jsn.cap);
            sn.classcap = JLINE.from_jline_matrix(jsn.classcap);
            sn.refclass = JLINE.from_jline_matrix(jsn.refclass)+1;
            sn.lldscaling = JLINE.from_jline_matrix(jsn.lldscaling);

            if ~isempty(jsn.cdscaling) && jsn.cdscaling.size() > 0
                % Convert Java SerializableFunction to MATLAB function handles
                sn.cdscaling = cell(sn.nstations, 1);
                % Iterate through the map entries to handle null values properly
                entrySet = jsn.cdscaling.entrySet();
                entryIter = entrySet.iterator();
                stationFunMap = containers.Map();
                while entryIter.hasNext()
                    entry = entryIter.next();
                    stationName = char(entry.getKey().getName());
                    try
                        jfun = entry.getValue();
                        if ~isempty(jfun)
                            stationFunMap(stationName) = jfun;
                        end
                    catch
                        % getValue() returns null for default lambda functions
                        % Skip and use default value
                    end
                end
                % Assign functions to stations
                for i = 1:sn.nstations
                    jstation = jstations.get(i-1);
                    stationName = char(jstation.getName());
                    if isKey(stationFunMap, stationName)
                        jfun = stationFunMap(stationName);
                        % Create a MATLAB function handle that calls the Java apply() method
                        sn.cdscaling{i} = @(ni) JLINE.call_java_cdscaling(jfun, ni);
                    else
                        sn.cdscaling{i} = @(ni) 1;
                    end
                end
            else
                sn.cdscaling = cell(sn.nstations, 0);
            end

            if ~isempty(jsn.nodetype)
                sn.nodetype = zeros(sn.nnodes, 1);
                for i = 1:jsn.nodetype.size
                    nodetype = jsn.nodetype.get(i-1);
                    switch nodetype.name().toCharArray'
                        case 'Queue'
                            sn.nodetype(i) = NodeType.Queue;
                        case 'Delay'
                            sn.nodetype(i) = NodeType.Delay;
                        case 'Source'
                            sn.nodetype(i) = NodeType.Source;
                        case 'Sink'
                            sn.nodetype(i) = NodeType.Sink;
                        case 'Join'
                            sn.nodetype(i) = NodeType.Join;
                        case 'Fork'
                            sn.nodetype(i) = NodeType.Fork;
                        case 'ClassSwitch'
                            sn.nodetype(i) = NodeType.ClassSwitch;
                        case 'Logger'
                            sn.nodetype(i) = NodeType.Logger;
                        case 'Cache'
                            sn.nodetype(i) = NodeType.Cache;
                        case 'Place'
                            sn.nodetype(i) = NodeType.Place;
                        case 'Transition'
                            sn.nodetype(i) = NodeType.Transition;
                        case 'Router'
                            sn.nodetype(i) = NodeType.Router;
                    end
                end
            else
                sn.nodetype = [];
            end

            if ~isempty(jsn.classnames)
                for i = 1:jsn.classnames.size
                    sn.classnames(i,1) = jsn.classnames.get(i-1);
                end
            else
                sn.classnames = [];
            end

            if ~isempty(jsn.nodenames)
                for i = 1:jsn.nodenames.size
                    sn.nodenames(i,1) = jsn.nodenames.get(i-1);
                end
            else
                sn.nodenames = [];
            end

            if ~isempty(jsn.rtorig) && jsn.rtorig.size()>0
                sn.rtorig = cell(sn.nclasses, sn.nclasses);
                for r = 1:sn.nclasses
                    for s = 1:sn.nclasses
                        sn.rtorig{r,s} = JLINE.from_jline_matrix(jsn.rtorig.get(jclasses.get(r-1)).get(jclasses.get(s-1)));
                    end
                end
            else
                sn.rtorig = {};
            end

            if ~isempty(jsn.state)
                sn.state = cell(sn.nstateful, 1);
                for i = 1:sn.nstateful
                    sn.state{i} = JLINE.from_jline_matrix(jstateful.get(i-1).getState());
                end
            else
                sn.state = {};
            end

            if ~isempty(jsn.stateprior)
                sn.stateprior = cell(sn.nstateful, 1);
                for i = 1:sn.nstateful
                    sn.stateprior{i} = JLINE.from_jline_matrix(jstateful.get(i-1).getStatePrior());
                end
            else
                sn.stateprior = {};
            end

            if ~isempty(jsn.space)
                sn.space = cell(sn.nstateful, 1);
                for i = 1:sn.nstateful
                    sn.space{i} = JLINE.from_jline_matrix(jstateful.get(i-1).getStateSpace());
                end
            else
                sn.space = {};
            end

            if ~isempty(jsn.routing)
                sn.routing = zeros(sn.nnodes, sn.nclasses);
                for i = 1:sn.nnodes
                    for j = 1:sn.nclasses
                        routingStrategy = jsn.routing.get(jnodes.get(i-1)).get(jclasses.get(j-1));
                        switch routingStrategy.name().toCharArray'
                            case 'PROB'
                                sn.routing(i,j) = RoutingStrategy.PROB;
                            case 'RAND'
                                sn.routing(i,j) = RoutingStrategy.RAND;
                            case 'RROBIN'
                                sn.routing(i,j) = RoutingStrategy.RROBIN;
                            case 'WRROBIN'
                                sn.routing(i,j) = RoutingStrategy.WRROBIN;
                            case 'JSQ'
                                sn.routing(i,j) = RoutingStrategy.JSQ;
                            case 'DISABLED'
                                sn.routing(i,j) = RoutingStrategy.DISABLED;
                            case 'FIRING'
                                sn.routing(i,j) = RoutingStrategy.FIRING;
                            case 'KCHOICES'
                                sn.routing(i,j) = RoutingStrategy.KCHOICES;
                        end
                    end
                end
            else
                sn.routing = [];
            end

            if ~isempty(jsn.procid)
                sn.procid = nan(sn.nstations, sn.nclasses);  % Initialize with NaN to match MATLAB behavior
                for i = 1:sn.nstations
                    for j = 1:sn.nclasses
                        stationMap = jsn.procid.get(jstations.get(i-1));
                        if isempty(stationMap)
                            sn.procid(i,j) = ProcessType.DISABLED;
                            continue;
                        end
                        processType = stationMap.get(jclasses.get(j-1));
                        if isempty(processType)
                            sn.procid(i,j) = ProcessType.DISABLED;
                            continue;
                        end
                        switch processType.name.toCharArray'
                            case 'EXP'
                                sn.procid(i,j) = ProcessType.EXP;
                            case 'ERLANG'
                                sn.procid(i,j) = ProcessType.ERLANG;
                            case 'HYPEREXP'
                                sn.procid(i,j) = ProcessType.HYPEREXP;
                            case 'PH'
                                sn.procid(i,j) = ProcessType.PH;
                            case 'APH'
                                sn.procid(i,j) = ProcessType.APH;
                            case 'MAP'
                                sn.procid(i,j) = ProcessType.MAP;
                            case 'UNIFORM'
                                sn.procid(i,j) = ProcessType.UNIFORM;
                            case 'DET'
                                sn.procid(i,j) = ProcessType.DET;
                            case 'COXIAN'
                                sn.procid(i,j) = ProcessType.COXIAN;
                            case 'GAMMA'
                                sn.procid(i,j) = ProcessType.GAMMA;
                            case 'PARETO'
                                sn.procid(i,j) = ProcessType.PARETO;
                            case 'WEIBULL'
                                sn.procid(i,j) = ProcessType.WEIBULL;
                            case 'LOGNORMAL'
                                sn.procid(i,j) = ProcessType.LOGNORMAL;
                            case 'MMPP2'
                                sn.procid(i,j) = ProcessType.MMPP2;
                            case 'REPLAYER'
                                sn.procid(i,j) = ProcessType.REPLAYER;
                            case 'TRACE'
                                sn.procid(i,j) = ProcessType.TRACE;
                            case 'IMMEDIATE'
                                sn.procid(i,j) = ProcessType.IMMEDIATE;
                            case 'DISABLED'
                                sn.procid(i,j) = ProcessType.DISABLED;
                            case 'COX2'
                                sn.procid(i,j) = ProcessType.COX2;
                            case 'BMAP'
                                sn.procid(i,j) = ProcessType.BMAP;
                            case 'ME'
                                sn.procid(i,j) = ProcessType.ME;
                            case 'RAP'
                                sn.procid(i,j) = ProcessType.RAP;
                            case 'BINOMIAL'
                                sn.procid(i,j) = ProcessType.BINOMIAL;
                            case 'POISSON'
                                sn.procid(i,j) = ProcessType.POISSON;
                            case 'GEOMETRIC'
                                sn.procid(i,j) = ProcessType.GEOMETRIC;
                            case 'DUNIFORM'
                                sn.procid(i,j) = ProcessType.DUNIFORM;
                            case 'BERNOULLI'
                                sn.procid(i,j) = ProcessType.BERNOULLI;
                            case 'PRIOR'
                                sn.procid(i,j) = ProcessType.PRIOR;
                            otherwise
                                % Unknown ProcessType - default to DISABLED
                                sn.procid(i,j) = ProcessType.DISABLED;
                        end
                    end
                end
            else
                sn.procid = [];
            end

            if ~isempty(jsn.mu)
                sn.mu = cell(sn.nstations, 1);
                for i = 1:sn.nstations
                    sn.mu{i} = cell(1, sn.nclasses);
                    for j = 1:sn.nclasses
                        sn.mu{i}{j} = JLINE.from_jline_matrix(jsn.mu.get(jstations.get(i-1)).get(jclasses.get(j-1)));
                    end
                end
            else
                sn.mu = {};
            end

            if ~isempty(jsn.phi)
                sn.phi = cell(sn.nstations, 1);
                for i = 1:sn.nstations
                    sn.phi{i} = cell(1, sn.nclasses);
                    for j = 1:sn.nclasses
                        sn.phi{i}{j} = JLINE.from_jline_matrix(jsn.phi.get(jstations.get(i-1)).get(jclasses.get(j-1)));
                    end
                end
            else
                sn.phi = {};
            end

            if ~isempty(jsn.proc)
                sn.proc = cell(sn.nstations, 1);
                for i = 1:sn.nstations
                    sn.proc{i} = cell(1, sn.nclasses);
                    for j = 1:sn.nclasses
                        proc_i_j = jsn.proc.get(jstations.get(i-1)).get(jclasses.get(j-1));
                        sn.proc{i}{j} = cell(1, proc_i_j.size);
                        for k = 1:proc_i_j.size
                            sn.proc{i}{j}{k} = JLINE.from_jline_matrix(proc_i_j.get(uint32(k-1)));
                        end
                    end
                end
            else
                sn.proc = {};
            end

            if ~isempty(jsn.pie)
                sn.pie = cell(sn.nstations, 1);
                for i = 1:sn.nstations
                    sn.pie{i} = cell(1, sn.nclasses);
                    for j = 1:sn.nclasses
                        sn.pie{i}{j} = JLINE.from_jline_matrix(jsn.pie.get(jstations.get(i-1)).get(jclasses.get(j-1)));
                    end
                end
            else
                sn.pie = {};
            end

            if ~isempty(jsn.sched)
                sn.sched = zeros(sn.nstations, 1);
                for i = 1:sn.nstations
                    schedStrategy = jsn.sched.get(jstations.get(i-1));
                    switch schedStrategy.name.toCharArray'
                        case 'INF'
                            sn.sched(i) = SchedStrategy.INF;
                        case 'FCFS'
                            sn.sched(i) = SchedStrategy.FCFS;
                        case 'LCFS'
                            sn.sched(i) = SchedStrategy.LCFS;
                        case 'LCFSPR'
                            sn.sched(i) = SchedStrategy.LCFSPR;
                        case 'SIRO'
                            sn.sched(i) = SchedStrategy.SIRO;
                        case 'SJF'
                            sn.sched(i) = SchedStrategy.SJF;
                        case 'LJF'
                            sn.sched(i) = SchedStrategy.LJF;
                        case 'PS'
                            sn.sched(i) = SchedStrategy.PS;
                        case 'DPS'
                            sn.sched(i) = SchedStrategy.DPS;
                        case 'GPS'
                            sn.sched(i) = SchedStrategy.GPS;
                        case 'PSPRIO'
                            sn.sched(i) = SchedStrategy.PSPRIO;
                        case 'DPSPRIO'
                            sn.sched(i) = SchedStrategy.DPSPRIO;
                        case 'GPSPRIO'
                            sn.sched(i) = SchedStrategy.GPSPRIO;
                        case 'SEPT'
                            sn.sched(i) = SchedStrategy.SEPT;
                        case 'LEPT'
                            sn.sched(i) = SchedStrategy.LEPT;
                        case {'HOL', 'FCFSPRIO'}
                            sn.sched(i) = SchedStrategy.FCFSPRIO;
                        case 'FORK'
                            sn.sched(i) = SchedStrategy.FORK;
                        case 'EXT'
                            sn.sched(i) = SchedStrategy.EXT;
                        case 'REF'
                            sn.sched(i) = SchedStrategy.REF;
                    end
                end
            else
                sn.sched = [];
            end

            if ~isempty(jsn.inchain)
                sn.inchain = cell(1, sn.nchains);
                for i = 1:sn.nchains
                    sn.inchain{1,i} = JLINE.from_jline_matrix(jsn.inchain.get(uint32(i-1)))+1;
                end
            else
                sn.inchain = {};
            end

            if ~isempty(jsn.visits)
                sn.visits = cell(sn.nchains, 1);
                for i = 1:sn.nchains
                    sn.visits{i,1} = JLINE.from_jline_matrix(jsn.visits.get(uint32(i-1)));
                end
            else
                sn.visits = {};
            end

            if ~isempty(jsn.nodevisits)
                sn.nodevisits = cell(1, sn.nchains);
                for i = 1:sn.nchains
                    sn.nodevisits{1,i} = JLINE.from_jline_matrix(jsn.nodevisits.get(uint32(i-1)));
                end
            else
                sn.nodevisits = {};
            end

            if ~isempty(jsn.droprule)
                sn.droprule = zeros(sn.nstations, sn.nclasses);
                for i = 1:sn.nstations
                    for j = 1:sn.nclasses
                        dropStrategy = jsn.droprule.get(jstations.get(i-1)).get(jclasses.get(j-1));
                        switch dropStrategy.name.toCharArray'
                            case 'WaitingQueue'
                                sn.droprule(i,j) = DropStrategy.WAITQ;
                            case 'Drop'
                                sn.droprule(i,j) = DropStrategy.DROP;
                            case 'BlockingAfterService'
                                sn.droprule(i,j) = DropStrategy.BAS;
                        end
                    end
                end
            else
                sn.droprule = [];
            end

            if ~isempty(jsn.nodeparam)
                sn.nodeparam = cell(sn.nnodes, 1);

                for i = 1:sn.nnodes
                    jnode = jnodes.get(i-1);
                    jparam = jsn.nodeparam.get(jnode);

                    %if jparam.isEmpty
                    %    sn.nodeparam{i} = [];
                    %    continue;
                    %end

                    % StationNodeParam
                    if isa(jparam, 'jline.lang.nodeparam.StationNodeParam')
                        if ~isempty(jparam.fileName)
                            sn.nodeparam{i}.fileName = cell(1, sn.nclasses);
                            for r = 1:sn.nclasses
                                fname = jparam.fileName.get(r-1);
                                if ~isempty(fname)
                                    sn.nodeparam{i}.fileName{r} = char(fname);
                                end
                            end
                        end
                    end

                    % TransitionNodeParam
                    if isa(jparam, 'jline.lang.nodeparam.TransitionNodeParam')
                        if ~isempty(jparam.firingprocid)
                            sn.nodeparam{i}.firingprocid = containers.Map('KeyType', 'char', 'ValueType', 'any');
                            keys = jparam.firingprocid.keySet.iterator;
                            while keys.hasNext
                                key = keys.next;
                                proc = jparam.firingprocid.get(key);
                                sn.nodeparam{i}.firingprocid(char(key.toString)) = char(proc.toString);
                            end
                        end
                        if ~isempty(jparam.firingphases)
                            sn.nodeparam{i}.firingphases = JLINE.from_jline_matrix(jparam.firingphases);
                        end
                        if ~isempty(jparam.fireweight)
                            sn.nodeparam{i}.fireweight = JLINE.from_jline_matrix(jparam.fireweight);
                        end
                    end

                    % JoinNodeParam
                    if isa(jparam, 'jline.lang.nodeparam.JoinNodeParam')
                        if ~isempty(jparam.joinStrategy)
                            sn.nodeparam{i}.joinStrategy = cell(1, sn.nclasses);
                            sn.nodeparam{i}.fanIn = cell(1, sn.nclasses);
                            for r = 1:sn.nclasses
                                jclass = jclasses.get(r-1);
                                joinStrategy = jparam.joinStrategy.get(jclass);
                                if ~isempty(joinStrategy)
                                    strategyStr = char(joinStrategy.name.toString);
                                    switch strategyStr
                                        case 'STD'
                                            sn.nodeparam{i}.joinStrategy{r} = JoinStrategy.STD;
                                        case 'PARTIAL'
                                            sn.nodeparam{i}.joinStrategy{r} = JoinStrategy.PARTIAL;
                                        otherwise
                                            sn.nodeparam{i}.joinStrategy{r} = strategyStr;
                                    end
                                    sn.nodeparam{i}.fanIn{r} = jparam.fanIn.get(jclass);
                                end
                            end
                        end
                    end

                    % RoutingNodeParam
                    if isa(jparam, 'jline.lang.nodeparam.RoutingNodeParam')
                        for r = 1:sn.nclasses
                            jclass = jclasses.get(r-1);

                            if ~isempty(jparam.weights) && jparam.weights.containsKey(jclass)
                                sn.nodeparam{i}.weights{r} = JLINE.from_jline_matrix(jparam.weights.get(jclass));
                            end

                            if ~isempty(jparam.outlinks) && jparam.outlinks.containsKey(jclass)
                                sn.nodeparam{i}.outlinks{r} = JLINE.from_jline_matrix(jparam.outlinks.get(jclass));
                            end
                        end
                    end

                    % ForkNodeParam
                    if isa(jparam, 'jline.lang.nodeparam.ForkNodeParam')
                        if ~isnan(jparam.fanOut)
                            sn.nodeparam{i}.fanOut = jparam.fanOut;
                        end
                    end

                    % CacheNodeParam
                    if isa(jparam, 'jline.lang.nodeparam.CacheNodeParam')
                        % nitems
                        if ~isnan(jparam.nitems)
                            sn.nodeparam{i}.nitems = jparam.nitems;
                        end

                        % accost
                        if ~isempty(jparam.accost)
                            % For Java 2D arrays (Matrix[][]), size(arr,2) returns 1 in MATLAB
                            % We need to get length of first row to get actual second dimension
                            K1 = size(jparam.accost, 1);
                            if K1 > 0
                                firstRow = jparam.accost(1);  % Get first row (Java array)
                                K2 = length(firstRow);
                            else
                                K2 = 0;
                            end
                            sn.nodeparam{i}.accost = cell(K1, K2);
                            for k1 = 1:K1
                                for k2 = 1:K2
                                    mat = jparam.accost(k1, k2); % MATLAB handles Java array indexing
                                    if ~isempty(mat)
                                        sn.nodeparam{i}.accost{k1, k2} = JLINE.from_jline_matrix(mat);
                                    end
                                end
                            end
                        end

                        % itemcap
                        if ~isempty(jparam.itemcap)
                            sn.nodeparam{i}.itemcap = JLINE.from_jline_matrix(jparam.itemcap);
                        end

                        % pread - convert from Java Map<Integer, List<Double>> to MATLAB cell array {R}
                        if ~isempty(jparam.pread)
                            nclasses = sn.nclasses;
                            sn.nodeparam{i}.pread = cell(1, nclasses);
                            for r = 1:nclasses
                                list = jparam.pread.get(int32(r-1)); % Java 0-based indexing
                                if ~isempty(list)
                                    values = zeros(1, list.size);
                                    for j = 1:list.size
                                        values(j) = list.get(j-1);
                                    end
                                    sn.nodeparam{i}.pread{r} = values;
                                else
                                    sn.nodeparam{i}.pread{r} = NaN;
                                end
                            end
                        end

                        % replacestrat
                        if ~isempty(jparam.replacestrat)
                            switch char(jparam.replacestrat)
                                case 'RR'
                                    sn.nodeparam{i}.replacestrat = ReplacementStrategy.RR;
                                case 'FIFO'
                                    sn.nodeparam{i}.replacestrat = ReplacementStrategy.FIFO;
                                case 'SFIFO'
                                    sn.nodeparam{i}.replacestrat = ReplacementStrategy.SFIFO;
                                case 'LRU'
                                    sn.nodeparam{i}.replacestrat = ReplacementStrategy.LRU;
                            end
                        end

                        % hitclass
                        if ~isempty(jparam.hitclass)
                            sn.nodeparam{i}.hitclass = 1+JLINE.from_jline_matrix(jparam.hitclass);
                        end
                        
                        % missclass
                        if ~isempty(jparam.missclass)
                            sn.nodeparam{i}.missclass =1+ JLINE.from_jline_matrix(jparam.missclass);
                        end

                        % actual hit/miss probabilities
                        if ~isempty(jparam.actualhitprob)
                            sn.nodeparam{i}.actualhitprob = JLINE.from_jline_matrix(jparam.actualhitprob);
                        end
                        if ~isempty(jparam.actualmissprob)
                            sn.nodeparam{i}.actualmissprob = JLINE.from_jline_matrix(jparam.actualmissprob);
                        end
                    end
                end
            else
                sn.nodeparam = {};
            end

            % if ~isempty(jsn.nodeparam)
            %     sn.nodeparam = cell(sn.nnodes, 1);
            %     % Note that JLINE only support node parameters related to
            %     % Fork, Join, WWROBIN and RROBIN
            %     for i = 1:sn.nnodes
            %         if jsn.nodeparam.get(jnodes.get(i-1)).isEmpty
            %             sn.nodeparam{i} = [];
            %         else
            %             if ~isnan(jsn.nodeparam.get(jnodes.get(i-1)).nitems)
            %                 sn.nodeparam{i}.nitems = jsn.nodeparam.get(jnodes.get(i-1)).nitems;
            %             end
            %             if ~isnan(jsn.nodeparam.get(jnodes.get(i-1)).fanOut)
            %                 sn.nodeparam{i}.fanOut = jsn.nodeparam.get(jnodes.get(i-1)).fanOut;
            %             end
            %             if ~isempty(jsn.nodeparam.get(jnodes.get(i-1)).joinStrategy)
            %                 if ~jsn.nodeparam.get(jnodes.get(i-1)).joinStrategy.isEmpty
            %                     sn.nodeparam{i}.joinStrategy = cell(1, sn.nclasses);
            %                     sn.nodeparam{i}.fanIn = cell(1, sn.nclasses);
            %                     for r = 1:sn.nclasses
            %                         joinStrategy = jsn.nodeparam.get(jnodes.get(i-1)).joinStrategy.get(jclasses.get(r-1));
            %                         switch joinStrategy.name.toCharArray'
            %                             case 'STD'
            %                                 sn.nodeparam{i}.joinStrategy{r} = JoinStrategy.STD;
            %                             case 'PARTIAL'
            %                                 sn.nodeparam{i}.joinStrategy{r} = JoinStrategy.PARTIAL;
            %                         end
            %                         sn.nodeparam{i}.fanIn{r} = jsn.nodeparam.get(jnodes.get(i-1)).fanIn.get(jclasses.get(r-1));
            %                     end
            %                 end
            %             end
            %
            %             if ~isempty(jsn.nodeparam.get(jnodes.get(i-1)).weights)
            %                 for r = 1:sn.nclasses
            %                     sn.nodeparam{i}{r}.weights = JLINE.from_jline_matrix(jsn.nodeparam.get(jnodes.get(i-1)).weights.get(jclasses.get(r-1)));
            %                 end
            %             end
            %
            %             if ~isempty(jsn.nodeparam.get(jnodes.get(i-1)).outlinks)
            %                 for r = 1:sn.nclasses
            %                     sn.nodeparam{i}{r}.outlinks = JLINE.from_jline_matrix(jsn.nodeparam.get(jnodes.get(i-1)).outlinks.get(jclasses.get(r-1)));
            %                 end
            %             end
            %         end
            %     end
            % else
            %     sn.nodeparam = {};
            % end

            if ~isempty(jsn.sync)
                jsync = jsn.sync;
                sn.sync = cell(jsync.size, 1);
                for i = 1:jsync.size
                    jsync_i = jsync.get(uint32(i-1));
                    sn.sync{i,1} = struct('active',cell(1),'passive',cell(1));

                    jactive = jsync_i.active.get(uint32(0));
                    jpassive = jsync_i.passive.get(uint32(0));

                    %Currently assume that prob would always be a value
                    %instead of lambda function (No idea of how to convert
                    %Java lambda function to matlab lambda function)
                    switch jactive.getEvent.name.toCharArray'
                        case 'INIT'
                            sn.sync{i,1}.active{1} = Event(EventType.INIT, jactive.getNode+1, jactive.getJobClass+1, ...
                                jactive.getProb, JLINE.from_jline_matrix(jactive.getState), ...
                                jactive.getT, jactive.getJob);
                        case 'LOCAL'
                            sn.sync{i,1}.active{1} = Event(EventType.LOCAL, jactive.getNode+1, jactive.getJobClass+1, ...
                                jactive.getProb, JLINE.from_jline_matrix(jactive.getState), ...
                                jactive.getT, jactive.getJob);
                        case 'ARV'
                            sn.sync{i,1}.active{1} = Event(EventType.ARV, jactive.getNode+1, jactive.getJobClass+1, ...
                                jactive.getProb, JLINE.from_jline_matrix(jactive.getState), ...
                                jactive.getT, jactive.getJob);
                        case 'DEP'
                            sn.sync{i,1}.active{1} = Event(EventType.DEP, jactive.getNode+1, jactive.getJobClass+1, ...
                                jactive.getProb, JLINE.from_jline_matrix(jactive.getState), ...
                                jactive.getT, jactive.getJob);
                        case 'PHASE'
                            sn.sync{i,1}.active{1} = Event(EventType.PHASE, jactive.getNode+1, jactive.getJobClass+1, ...
                                jactive.getProb, JLINE.from_jline_matrix(jactive.getState), ...
                                jactive.getT, jactive.getJob);
                        case 'READ'
                            sn.sync{i,1}.active{1} = Event(EventType.READ, jactive.getNode+1, jactive.getJobClass+1, ...
                                jactive.getProb, JLINE.from_jline_matrix(jactive.getState), ...
                                jactive.getT, jactive.getJob);
                        case 'STAGE'
                            sn.sync{i,1}.active{1} = Event(EventType.STAGE, jactive.getNode+1, jactive.getJobClass+1, ...
                                jactive.getProb, JLINE.from_jline_matrix(jactive.getState), ...
                                jactive.getT, jactive.getJob);
                    end

                    switch jpassive.getEvent.name.toCharArray'
                        case 'INIT'
                            sn.sync{i,1}.passive{1} = Event(EventType.INIT, jpassive.getNode+1, jpassive.getJobClass+1, ...
                                jpassive.getProb, JLINE.from_jline_matrix(jpassive.getState), ...
                                jpassive.getT, jpassive.getJob);
                        case 'LOCAL'
                            sn.sync{i,1}.passive{1} = Event(EventType.LOCAL, jpassive.getNode+1, jpassive.getJobClass+1, ...
                                jpassive.getProb, JLINE.from_jline_matrix(jpassive.getState), ...
                                jpassive.getT, jpassive.getJob);
                        case 'ARV'
                            sn.sync{i,1}.passive{1} = Event(EventType.ARV, jpassive.getNode+1, jpassive.getJobClass+1, ...
                                jpassive.getProb, JLINE.from_jline_matrix(jpassive.getState), ...
                                jpassive.getT, jpassive.getJob);
                        case 'DEP'
                            sn.sync{i,1}.passive{1} = Event(EventType.DEP, jpassive.getNode+1, jpassive.getJobClass+1, ...
                                jpassive.getProb, JLINE.from_jline_matrix(jpassive.getState), ...
                                jpassive.getT, jpassive.getJob);
                        case 'PHASE'
                            sn.sync{i,1}.passive{1} = Event(EventType.PHASE, jpassive.getNode+1, jpassive.getJobClass+1, ...
                                jpassive.getProb, JLINE.from_jline_matrix(jpassive.getState), ...
                                jpassive.getT, jpassive.getJob);
                        case 'READ'
                            sn.sync{i,1}.passive{1} = Event(EventType.READ, jpassive.getNode+1, jpassive.getJobClass+1, ...
                                jpassive.getProb, JLINE.from_jline_matrix(jpassive.getState), ...
                                jpassive.getT, jpassive.getJob);
                        case 'STAGE'
                            sn.sync{i,1}.passive{1} = Event(EventType.STAGE, jpassive.getNode+1, jpassive.getJobClass+1, ...
                                jpassive.getProb, JLINE.from_jline_matrix(jpassive.getState), ...
                                jpassive.getT, jpassive.getJob);
                    end
                end
            else
                sn.sync = {};
            end
        end

        function [QN,UN,RN,WN,AN,TN] = arrayListToResults(alist)
            switch class(alist)
                case 'jline.solvers.LayeredNetworkAvgTable'
                    QN = JLINE.arraylist_to_matrix(alist.getQLen());
                    UN = JLINE.arraylist_to_matrix(alist.getUtil());
                    RN = JLINE.arraylist_to_matrix(alist.getRespT());
                    WN = JLINE.arraylist_to_matrix(alist.getResidT());
                    AN = NaN*JLINE.arraylist_to_matrix(alist.getTput()); % getArvR not yet available in JLINE
                    TN = JLINE.arraylist_to_matrix(alist.getTput());
                otherwise
                    QN = JLINE.arraylist_to_matrix(alist.getQLen());
                    UN = JLINE.arraylist_to_matrix(alist.getUtil());
                    RN = JLINE.arraylist_to_matrix(alist.getRespT());
                    WN = JLINE.arraylist_to_matrix(alist.getResidT());
                    AN = JLINE.arraylist_to_matrix(alist.getArvR());
                    TN = JLINE.arraylist_to_matrix(alist.getTput());
            end
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source',...
                'ClassSwitch','Delay','DelayStation','Queue',...
                'APH','Coxian','Erlang','Exp','HyperExp',...
                'StatelessClassSwitcher','InfiniteServer','SharedServer','Buffer','Dispatcher',...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ClosedClass','OpenClass'});
        end

        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)

            featUsed = model.getUsedLangFeatures();
            featSupported = JLINE.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end


        function solverOptions = parseSolverOptions(solverOptions, options)
            fn = fieldnames(options);
            fn2 = fieldnames(solverOptions);
            for f = 1:length(fn)
                found = 0;
                for j = 1:length(fn2)
                    if strcmp(fn{f}, fn2{j})
                        found = 1;
                        switch fn{f}
                            case 'seed'
                                solverOptions.seed = options.seed;
                            case 'samples'
                                solverOptions.samples = options.samples;
                            case 'confint'
                                % Parse confint - can be a level (0.95) or 0 to disable
                                [confintEnabled, confintLevel] = Solver.parseConfInt(options.confint);
                                if confintEnabled
                                    solverOptions.confint = confintLevel;
                                else
                                    solverOptions.confint = 0;
                                end
                            case 'method'
                                solverOptions.method = options.method;
                            case 'config' % SSA specific
                                if isfield(options.config,'eventcache')
                                    solverOptions.config.eventcache = options.config.eventcache;
                                end
                            case 'verbose'
                                switch options.(fn{f})
                                    case {VerboseLevel.SILENT}
                                        solverOptions.verbose = solverOptions.verbose.SILENT;
                                    case {VerboseLevel.STD}
                                        solverOptions.verbose = solverOptions.verbose.STD;
                                    case {VerboseLevel.DEBUG}
                                        solverOptions.verbose = solverOptions.verbose.DEBUG;
                                end
                            case 'init_sol'
                                solverOptions.(fn{f}) = JLINE.from_line_matrix(options.init_sol);
                            case 'cutoff'
                                if isscalar(options.cutoff)
                                    solverOptions.(fn{f}) = jline.util.matrix.Matrix.singleton(options.cutoff);
                                else
                                    solverOptions.(fn{f}) = JLINE.from_line_matrix(options.cutoff);
                                end
                            case 'odesolvers'
                            case 'rewardIterations'
                                solverOptions.rewardIterations = java.lang.Integer(options.rewardIterations);
                            otherwise
                                solverOptions.(fn{f}) = options.(fn{f});
                        end

                        break;
                    end
                end
                if ~found
                    line_printf('Could not find option %s in the JLINE options.\n', fn{f});
                end
            end
        end

        function [ssa] = SolverSSA(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.SSA);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            jline.util.Maths.setRandomNumbersMatlab(true);
            ssa = jline.solvers.ssa.SolverSSA(network_object, solverOptions);
        end

        function [mam] = SolverMAM(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.MAM);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            mam = jline.solvers.mam.SolverMAM(network_object, solverOptions);
        end

        function [jmt] = SolverJMT(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.JMT);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            jmt = jline.solvers.jmt.SolverJMT(network_object, solverOptions);
        end

        function [ctmc] = SolverCTMC(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.CTMC);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            ctmc = jline.solvers.ctmc.SolverCTMC(network_object,solverOptions);
        end

        function [fluid] = SolverFluid(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.FLUID);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            fluid = jline.solvers.fluid.SolverFluid(network_object, solverOptions);
        end

        function [QN, UN, RN, TN, CN, XN, t, QNt, UNt, TNt, xvec] = runFluidAnalyzer(network, options)
            % RUNFLUIDANALYZER Run JLINE fluid analyzer and return results
            %
            % [QN, UN, RN, TN, CN, XN, T, QNT, UNT, TNT, XVEC] = JLINE.runFluidAnalyzer(NETWORK, OPTIONS)
            %
            % Runs the JLINE fluid solver on the given network and converts
            % results back to MATLAB data structures.
            %
            % Input:
            %   network - LINE Network model
            %   options - Solver options structure with fields:
            %             .method - solver method
            %             .stiff  - use stiff ODE solver
            %
            % Output:
            %   QN, UN, RN, TN - Steady-state metrics [M x K]
            %   CN, XN         - System metrics [1 x K]
            %   t              - Time vector [Tmax x 1]
            %   QNt, UNt, TNt  - Transient metrics {M x K} cells
            %   xvec           - State vector structure

            jmodel = LINE2JLINE(network);
            jsolver = JLINE.SolverFluid(jmodel);
            import jline.solvers.fluid.*;

            jsolver.options.method = options.method;
            jsolver.options.stiff = options.stiff;
            result = jsolver.runMethodSpecificAnalyzerViaLINE();

            % Convert JLINE result to MATLAB data structures
            M = jmodel.getNumberOfStatefulNodes();
            K = jmodel.getNumberOfClasses();

            QN = NaN * zeros(M, K);
            UN = NaN * zeros(M, K);
            RN = NaN * zeros(M, K);
            TN = NaN * zeros(M, K);
            CN = NaN * zeros(1, K);
            XN = NaN * zeros(1, K);

            QNt = cell(M, K);
            UNt = cell(M, K);
            TNt = cell(M, K);

            Tmax = result.t.length();
            t = NaN * zeros(Tmax, 1);

            for ist = 1:M
                for jst = 1:K
                    QN(ist, jst) = result.QN.get(ist-1, jst-1);
                    UN(ist, jst) = result.UN.get(ist-1, jst-1);
                    RN(ist, jst) = result.RN.get(ist-1, jst-1);
                    TN(ist, jst) = result.TN.get(ist-1, jst-1);
                end
            end

            for jst = 1:K
                CN(1, jst) = result.CN.get(0, jst-1);
                XN(1, jst) = result.XN.get(0, jst-1);
            end

            for ist = 1:M
                for jst = 1:K
                    for p = 1:Tmax
                        QNt{ist, jst}(p, 1) = result.QNt(ist, jst).get(p-1, 0);
                        UNt{ist, jst}(p, 1) = result.UNt(ist, jst).get(p-1, 0);
                        TNt{ist, jst}(p, 1) = result.TNt(ist, jst).get(p-1, 0);
                    end
                end
            end

            for p = 1:Tmax
                t(p, 1) = result.t.get(p-1, 0);
            end

            % JLINE does not return odeStateVec
            xvec.odeStateVec = [];
            xvec.sn = network;
        end

        function [des] = SolverDES(network_object, options)
            % Create DES-specific options object
            desOptions = jline.solvers.des.DESOptions();
            if nargin>1
                % Copy standard options
                desOptions.samples = options.samples;
                desOptions.seed = options.seed;
                % Parse confint
                [confintEnabled, confintLevel] = Solver.parseConfInt(options.confint);
                if confintEnabled
                    desOptions.confint = confintLevel;
                else
                    desOptions.confint = 0;
                end
                % Pass timespan for transient analysis
                if isfield(options, 'timespan') && length(options.timespan) >= 2
                    desOptions.timespan = options.timespan;
                end
                % Pass DES-specific options if configured
                if isfield(options, 'config')
                    % Transient detection options
                    if isfield(options.config, 'tranfilter')
                        desOptions.tranfilter = options.config.tranfilter;
                    end
                    if isfield(options.config, 'mserbatch')
                        desOptions.mserbatch = options.config.mserbatch;
                    end
                    if isfield(options.config, 'warmupfrac')
                        desOptions.warmupfrac = options.config.warmupfrac;
                    end
                    % Confidence interval options
                    if isfield(options.config, 'cimethod')
                        desOptions.cimethod = options.config.cimethod;
                    end
                    if isfield(options.config, 'obmoverlap')
                        desOptions.obmoverlap = options.config.obmoverlap;
                    end
                    if isfield(options.config, 'ciminbatch')
                        desOptions.ciminbatch = options.config.ciminbatch;
                    end
                    if isfield(options.config, 'ciminobs')
                        desOptions.ciminobs = options.config.ciminobs;
                    end
                    % Convergence options
                    if isfield(options.config, 'cnvgon')
                        desOptions.cnvgon = options.config.cnvgon;
                    end
                    if isfield(options.config, 'cnvgtol')
                        desOptions.cnvgtol = options.config.cnvgtol;
                    end
                    if isfield(options.config, 'cnvgbatch')
                        desOptions.cnvgbatch = options.config.cnvgbatch;
                    end
                    if isfield(options.config, 'cnvgchk')
                        desOptions.cnvgchk = options.config.cnvgchk;
                    end
                end
            end
            des = jline.solvers.des.SolverDES(network_object, desOptions);
        end

        function [mva] = SolverMVA(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.MVA);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            mva = jline.solvers.mva.SolverMVA(network_object, solverOptions);
        end

        function [nc] = SolverNC(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.NC);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            nc = jline.solvers.nc.SolverNC(network_object, solverOptions);
        end

        function [auto] = SolverAuto(network_object, options)
            solverOptions = jline.solvers.auto.AUTOptions();
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            auto = jline.solvers.auto.SolverAUTO(network_object, solverOptions);
        end

        function streamOpts = StreamingOptions(varargin)
            % STREAMINGOPTIONS Create Java StreamingOptions for SSA/DES stream() method
            %
            % @brief Creates StreamingOptions for streaming simulation metrics
            %
            % @param varargin Name-value pairs for options:
            %   'transport' - 'http' (recommended) or 'grpc' (default: 'http')
            %   'endpoint' - Receiver endpoint (default: 'localhost:8080/metrics' for HTTP)
            %   'mode' - 'sampled' or 'time_window' (default: 'sampled')
            %   'sampleFrequency' - Push every N events in sampled mode (default: 100)
            %   'timeWindowSeconds' - Window duration in time_window mode (default: 1.0)
            %   'serviceName' - Service identifier (default: 'line-stream')
            %   'includeQueueLength' - Include queue length metrics (default: true)
            %   'includeUtilization' - Include utilization metrics (default: true)
            %   'includeThroughput' - Include throughput metrics (default: true)
            %   'includeResponseTime' - Include response time metrics (default: true)
            %   'includeArrivalRate' - Include arrival rate metrics (default: true)
            %
            % @return streamOpts Java StreamingOptions object
            %
            % Example:
            % @code
            % streamOpts = JLINE.StreamingOptions('transport', 'http', 'sampleFrequency', 50);
            % @endcode

            streamOpts = jline.streaming.StreamingOptions();

            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'transport', 'http', @ischar);
            addParameter(p, 'endpoint', '', @ischar);  % Empty means use default for transport
            addParameter(p, 'mode', 'sampled', @ischar);
            addParameter(p, 'sampleFrequency', 100, @isnumeric);
            addParameter(p, 'timeWindowSeconds', 1.0, @isnumeric);
            addParameter(p, 'serviceName', 'line-stream', @ischar);
            addParameter(p, 'includeQueueLength', true, @islogical);
            addParameter(p, 'includeUtilization', true, @islogical);
            addParameter(p, 'includeThroughput', true, @islogical);
            addParameter(p, 'includeResponseTime', true, @islogical);
            addParameter(p, 'includeArrivalRate', true, @islogical);
            parse(p, varargin{:});

            % Set transport type
            transportTypes = javaMethod('values', 'jline.streaming.StreamingOptions$TransportType');
            switch lower(p.Results.transport)
                case 'http'
                    streamOpts.transport = transportTypes(1);  % HTTP
                case 'grpc'
                    streamOpts.transport = transportTypes(2);  % GRPC
                otherwise
                    streamOpts.transport = transportTypes(1);  % Default to HTTP
            end

            % Set endpoint (use provided or default based on transport)
            if ~isempty(p.Results.endpoint)
                streamOpts.endpoint = p.Results.endpoint;
            end
            % If empty, StreamingOptions uses its default for the transport type

            % Set mode
            streamModes = javaMethod('values', 'jline.streaming.StreamingOptions$StreamMode');
            switch lower(p.Results.mode)
                case 'sampled'
                    streamOpts.mode = streamModes(1);  % SAMPLED
                case 'time_window'
                    streamOpts.mode = streamModes(2);  % TIME_WINDOW
                otherwise
                    streamOpts.mode = streamModes(1);  % Default to SAMPLED
            end

            % Set other options
            streamOpts.sampleFrequency = p.Results.sampleFrequency;
            streamOpts.timeWindowSeconds = p.Results.timeWindowSeconds;
            streamOpts.serviceName = p.Results.serviceName;
            streamOpts.includeQueueLength = p.Results.includeQueueLength;
            streamOpts.includeUtilization = p.Results.includeUtilization;
            streamOpts.includeThroughput = p.Results.includeThroughput;
            streamOpts.includeResponseTime = p.Results.includeResponseTime;
            streamOpts.includeArrivalRate = p.Results.includeArrivalRate;
        end

        function result = convertSampleResult(jresult)
            % CONVERTSAMPLERESULT Convert Java sample result to MATLAB struct
            %
            % @brief Converts Java SampleNodeState to MATLAB structure
            %
            % @param jresult Java SampleNodeState object
            % @return result MATLAB struct with fields: t, state, isaggregate

            result = struct();

            % Convert time matrix
            if ~isempty(jresult.t)
                result.t = JLINE.from_jline_matrix(jresult.t);
            else
                result.t = [];
            end

            % Convert state matrix
            if ~isempty(jresult.state) && isa(jresult.state, 'jline.util.matrix.Matrix')
                result.state = JLINE.from_jline_matrix(jresult.state);
            else
                result.state = [];
            end

            result.isaggregate = jresult.isaggregate;
        end

        function [ln] = SolverLN(layered_network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.LN);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            ln = jline.solvers.ln.SolverLN(layered_network_object, solverOptions);
        end

        function serfun = handle_to_serializablefun(handle, sn)
            % HANDLE_TO_SERIALIZABLEFUN Convert MATLAB function handle to Java SerializableFunction
            %
            % This function pre-computes the function values for all possible state
            % combinations and creates a Java PrecomputedCDFunction object.
            %
            % @param handle MATLAB function handle that takes a vector ni and returns a scalar
            % @param sn Network struct containing njobs (population per class)
            % @return serfun Java PrecomputedCDFunction object

            % Get number of classes and maximum populations
            nclasses = sn.nclasses;
            njobs = sn.njobs;  % Population per class

            % For open classes (njobs=0), use a reasonable bound
            maxPop = njobs;
            for r = 1:nclasses
                if maxPop(r) == 0 || isinf(maxPop(r))
                    % For open classes, use sum of closed class populations or 100 as bound
                    maxPop(r) = max(100, sum(njobs(isfinite(njobs) & njobs > 0)));
                end
            end

            % Create Java PrecomputedCDFunction object
            serfun = jline.util.PrecomputedCDFunction(nclasses);

            % Enumerate all possible state combinations and pre-compute function values
            % Use recursive enumeration to handle arbitrary number of classes
            JLINE.enumerate_states(handle, serfun, maxPop, zeros(1, nclasses), 1);
        end

        function enumerate_states(handle, serfun, maxPop, currentState, classIdx)
            % ENUMERATE_STATES Recursively enumerate all state combinations
            %
            % @param handle MATLAB function handle
            % @param serfun Java PrecomputedCDFunction object to populate
            % @param maxPop Maximum population per class
            % @param currentState Current state being built
            % @param classIdx Current class index being enumerated

            nclasses = length(maxPop);

            if classIdx > nclasses
                % We have a complete state, compute and store the function value
                try
                    value = handle(currentState);
                    % Convert to Java int array and add to serfun
                    jstate = jline.util.matrix.Matrix(1, nclasses);
                    for r = 1:nclasses
                        jstate.set(0, r-1, currentState(r));
                    end
                    serfun.addValue(jstate, value);
                catch
                    % If function evaluation fails, skip this state
                end
                return;
            end

            % Enumerate all populations for current class
            for n = 0:maxPop(classIdx)
                currentState(classIdx) = n;
                JLINE.enumerate_states(handle, serfun, maxPop, currentState, classIdx + 1);
            end
        end

        function result = call_java_cdscaling(jfun, ni)
            % CALL_JAVA_CDSCALING Call a Java SerializableFunction for class dependence
            %
            % This function converts a MATLAB vector to a Java Matrix and calls
            % the Java function's apply() method.
            %
            % @param jfun Java SerializableFunction<Matrix, Double> object
            % @param ni MATLAB vector representing the state (jobs per class)
            % @return result The scaling factor returned by the Java function

            % Convert MATLAB vector to Java Matrix
            if isrow(ni)
                jmatrix = jline.util.matrix.Matrix(1, length(ni));
                for r = 1:length(ni)
                    jmatrix.set(0, r-1, ni(r));
                end
            else
                jmatrix = jline.util.matrix.Matrix(length(ni), 1);
                for r = 1:length(ni)
                    jmatrix.set(r-1, 0, ni(r));
                end
            end

            % Call the Java function and convert result to MATLAB double
            jresult = jfun.apply(jmatrix);
            result = double(jresult);
        end

    end
end
