classdef JLINE
    % JLINE Conversion utilities for JLINE format models
    %
    % JLINE provides static methods to convert between LINE MATLAB models 
    % and JLINE Java models. This class serves as the primary interface
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

    methods(Static)

        function jar_loc = get_jar_location()
            % Get jline.jar location, downloading if necessary.
            % Re-checks on each call, so deleted JAR triggers re-download.
            % Treats lang='java' as a wrapper that auto-downloads jline.jar if absent.
            jar_loc = which('jline.jar');
            if isempty(jar_loc) || ~isfile(jar_loc)
                jar_loc = lineDownloadJAR(false);  % Silent re-download
            end
        end

        function set_layered_rate_dependence(jelem, sn, eidx)
            % SET_LAYERED_RATE_DEPENDENCE(JELEM, SN, EIDX)
            %
            % Load dependence is a numeric vector and crosses to the JAR; the
            % class- and joint-dependent scalings are MATLAB function handles
            % with no Java counterpart, so they are refused rather than dropped.
            if isfield(sn,'lldscaling') && numel(sn.lldscaling)>=eidx && ~isempty(sn.lldscaling{eidx})
                alpha = sn.lldscaling{eidx};
                jalpha = javaObject('jline.util.matrix.Matrix', 1, numel(alpha));
                for jj = 1:numel(alpha)
                    jalpha.set(0, jj-1, alpha(jj));
                end
                jelem.setLoadDependence(jalpha);
            end
            hascd = isfield(sn,'cdscaling') && numel(sn.cdscaling)>=eidx && ~isempty(sn.cdscaling{eidx});
            hasjd = isfield(sn,'jdscaling') && numel(sn.jdscaling)>=eidx && ~isempty(sn.jdscaling{eidx});
            if hascd || hasjd
                line_error(mfilename,'Class- and joint-dependence on %s are MATLAB function handles and cannot be marshalled to the JAR; solve this model with lang=''matlab''.', sn.names{eidx});
            end
        end

        function set_layered_lincon(jelem, sn, eidx)
            % SET_LAYERED_LINCON(JELEM, SN, EIDX)
            %
            % Admission constraint A*n <= b of a Host or a Task. SN.LINCON holds
            % the resolved positional form, columns in the tasksof/entriesof order
            % that the JAR uses as well, so the named rows of ADDCONSTRAINT cross
            % already merged into A.
            if ~isfield(sn,'lincon') || size(sn.lincon,1) < eidx || isempty(sn.lincon{eidx,1})
                return
            end
            A = sn.lincon{eidx,1};
            b = sn.lincon{eidx,2};
            jA = javaObject('jline.util.matrix.Matrix', size(A,1), size(A,2));
            for ii = 1:size(A,1)
                for jj = 1:size(A,2)
                    jA.set(ii-1, jj-1, A(ii,jj));
                end
            end
            jb = javaObject('jline.util.matrix.Matrix', numel(b), 1);
            for ii = 1:numel(b)
                jb.set(ii-1, 0, b(ii));
            end
            jelem.setConstraint(jA, jb);
        end

        function model = from_line_layered_network(line_layered_network)
            sn = line_layered_network.getStruct;

            %% initialization
            model = javaObject('jline.lang.layered.LayeredNetwork', line_layered_network.getName);

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
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.REF);
                    case SchedStrategy.INF
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.INF);
                    case SchedStrategy.FCFS
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.FCFS);
                    case SchedStrategy.LCFS
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFS);
                    case SchedStrategy.SIRO
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.SIRO);
                    case SchedStrategy.SJF
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.SJF);
                    case SchedStrategy.LJF
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LJF);
                    case SchedStrategy.PS
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.PS);
                    case SchedStrategy.DPS
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.DPS);
                    case SchedStrategy.GPS
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.GPS);
                    case SchedStrategy.SEPT
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.SEPT);
                    case SchedStrategy.LEPT
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LEPT);
                    case {SchedStrategy.HOL, SchedStrategy.FCFSPRIO}
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.FCFSPRIO);
                    case SchedStrategy.FORK
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.FORK);
                    case SchedStrategy.EXT
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.EXT);
                    case SchedStrategy.LCFSPR
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFSPR);
                    case SchedStrategy.LCFSPI
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFSPI);
                    case SchedStrategy.LCFSPRIO
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFSPRIO);
                    case SchedStrategy.LCFSPRPRIO
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFSPRPRIO);
                    case SchedStrategy.LCFSPIPRIO
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.LCFSPIPRIO);
                    case SchedStrategy.PSPRIO % todo
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.PSPRIO);
                    case SchedStrategy.DPSPRIO % todo
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.DPSPRIO);
                    case SchedStrategy.GPSPRIO % todo
                        P{h} = javaObject('jline.lang.layered.Processor', model, sn.names{h}, sn_mult_h, jline.lang.constant.SchedStrategy.GPSPRIO);
                end
                if sn.repl(h)~=1
                    P{h}.setReplication(sn.repl(h));
                end
                JLINE.set_layered_rate_dependence(P{h}, sn, h);
                JLINE.set_layered_lincon(P{h}, sn, h);
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
                    T{t} = javaObject('jline.lang.layered.CacheTask', model, sn.names{tidx}, sn.nitems(tidx), sn.itemcap{tidx}, jReplacestrat, sn_mult_tidx);
                    if isfield(sn,'hasretrieval') && numel(sn.hasretrieval)>=tidx && sn.hasretrieval(tidx)
                        T{t}.setRetrieval(true);
                    end
                elseif sn.hassetup(tidx)
                    % SetupTask (has setupTime/delayOffTime)
                    jSchedStrategy = JLINE.to_jline_sched_strategy(sn.sched(tidx));
                    T{t} = javaObject('jline.lang.layered.SetupTask', model, sn.names{tidx}, sn_mult_tidx, jSchedStrategy);
                else
                    jSchedStrategy = JLINE.to_jline_sched_strategy(sn.sched(tidx));
                    T{t} = javaObject('jline.lang.layered.Task', model, sn.names{tidx}, sn_mult_tidx, jSchedStrategy);
                end
                T{t}.on(P{sn.parent(tidx)});
                if sn.repl(tidx)~=1
                    T{t}.setReplication(sn.repl(tidx));
                end
                JLINE.set_layered_rate_dependence(T{t}, sn, tidx);
                JLINE.set_layered_lincon(T{t}, sn, tidx);
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
                        case ProcessType.DET
                            T{t}.setThinkTime(jline.lang.processes.Det(sn.think_mean(tidx)));
                        case ProcessType.UNIFORM
                            p = sn.think_params{tidx};
                            T{t}.setThinkTime(jline.lang.processes.Uniform(p(1), p(2)));
                        case ProcessType.GAMMA
                            p = sn.think_params{tidx};
                            T{t}.setThinkTime(jline.lang.processes.Gamma(p(1), p(2)));
                        case ProcessType.PARETO
                            p = sn.think_params{tidx};
                            T{t}.setThinkTime(jline.lang.processes.Pareto(p(1), p(2)));
                        case ProcessType.WEIBULL
                            p = sn.think_params{tidx};
                            T{t}.setThinkTime(jline.lang.processes.Weibull(p(1), p(2)));
                        case ProcessType.LOGNORMAL
                            p = sn.think_params{tidx};
                            T{t}.setThinkTime(jline.lang.processes.Lognormal(p(1), p(2)));
                        otherwise
                            line_error(mfilename,sprintf('JLINE conversion does not support the %s distribution for task think time yet.',char(sn.think_type(tidx))));
                    end
                end
                % Setup time rebuilt from (type,mean,scv,params) -- see _kb/12-interfaces-and-docs.md
                if sn.hassetup(tidx) && ~isnan(sn.setuptime_mean(tidx)) && sn.setuptime_mean(tidx) > 1e-8
                    T{t}.setSetupTime(JLINE.from_line_lqn_dist(sn.setuptime_type(tidx), ...
                        sn.setuptime_mean(tidx), sn.setuptime_scv(tidx), ...
                        sn.setuptime_params{tidx}, sn.setuptime_proc{tidx}));
                end
                % Delay-off time
                if sn.hassetup(tidx) && ~isnan(sn.delayofftime_mean(tidx)) && sn.delayofftime_mean(tidx) > 1e-8
                    T{t}.setDelayOffTime(JLINE.from_line_lqn_dist(sn.delayofftime_type(tidx), ...
                        sn.delayofftime_mean(tidx), sn.delayofftime_scv(tidx), ...
                        sn.delayofftime_params{tidx}, sn.delayofftime_proc{tidx}));
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
                        jPopularity = javaObject('jline.lang.processes.DiscreteSampler', jline.util.matrix.Matrix.uniformDistribution(sn.nitems(eidx)));
                    end
                    E{e} = javaObject('jline.lang.layered.ItemEntry', model, sn.names{eidx}, sn.nitems(eidx), jPopularity);
                else
                    E{e} = javaObject('jline.lang.layered.Entry', model, sn.names{eidx});
                end
                E{e}.on(T{sn.parent(eidx)-sn.tshift});
                % Open arrival gated on arrival_mean, not arrival_type -- see _kb/12-interfaces-and-docs.md
                if ~isnan(sn.arrival_mean(eidx)) && sn.arrival_mean(eidx) > 0
                    E{e}.setArrival(JLINE.from_line_lqn_dist(sn.arrival_type(eidx), ...
                        sn.arrival_mean(eidx), sn.arrival_scv(eidx), ...
                        sn.arrival_params{eidx}, sn.arrival_proc{eidx}));
                end
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
                        jHostDem = javaObject('jline.lang.processes.Exp', 1/sn.hostdem_mean(aidx));
                    case ProcessType.ERLANG
                        jHostDem = jline.lang.processes.Erlang.fitMeanAndSCV(sn.hostdem_mean(aidx), sn.hostdem_scv(aidx));
                    case ProcessType.HYPEREXP
                        if ~isempty(sn.hostdem_params{aidx}) && length(sn.hostdem_params{aidx}) >= 3
                            jHostDem = javaObject('jline.lang.processes.HyperExp', sn.hostdem_params{aidx}(1), sn.hostdem_params{aidx}(2), sn.hostdem_params{aidx}(3));
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
                            jHostDem = javaObject('jline.lang.processes.PH', JLINE.from_line_matrix(proc{1}), JLINE.from_line_matrix(proc{2}));
                        else
                            jHostDem = javaObject('jline.lang.processes.Exp', 1/sn.hostdem_mean(aidx));
                        end
                    case ProcessType.MAP
                        if ~isempty(sn.hostdem_proc{aidx})
                            proc = sn.hostdem_proc{aidx};
                            jHostDem = javaObject('jline.lang.processes.MAP', JLINE.from_line_matrix(proc{1}), JLINE.from_line_matrix(proc{2}));
                        else
                            jHostDem = javaObject('jline.lang.processes.Exp', 1/sn.hostdem_mean(aidx));
                        end
                    case ProcessType.DET
                        jHostDem = javaObject('jline.lang.processes.Det', sn.hostdem_mean(aidx));
                    case ProcessType.UNIFORM
                        p = sn.hostdem_params{aidx};
                        jHostDem = javaObject('jline.lang.processes.Uniform', p(1), p(2));
                    case ProcessType.GAMMA
                        p = sn.hostdem_params{aidx};
                        jHostDem = javaObject('jline.lang.processes.Gamma', p(1), p(2));
                    case ProcessType.PARETO
                        p = sn.hostdem_params{aidx};
                        jHostDem = javaObject('jline.lang.processes.Pareto', p(1), p(2));
                    case ProcessType.WEIBULL
                        p = sn.hostdem_params{aidx};
                        jHostDem = javaObject('jline.lang.processes.Weibull', p(1), p(2));
                    case ProcessType.LOGNORMAL
                        p = sn.hostdem_params{aidx};
                        jHostDem = javaObject('jline.lang.processes.Lognormal', p(1), p(2));
                    otherwise
                        line_error(mfilename,sprintf('JLINE conversion does not support the %s distribution for host demand yet.',char(sn.hostdem_type(aidx))));
                end
                A{a} = javaObject('jline.lang.layered.Activity', model, sn.names{aidx}, jHostDem);
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
                        case ProcessType.DET
                            P{h}.setThinkTime(jline.lang.processes.Det(sn.think_mean(h)));
                        case ProcessType.UNIFORM
                            p = sn.think_params{h};
                            P{h}.setThinkTime(jline.lang.processes.Uniform(p(1), p(2)));
                        case ProcessType.GAMMA
                            p = sn.think_params{h};
                            P{h}.setThinkTime(jline.lang.processes.Gamma(p(1), p(2)));
                        case ProcessType.PARETO
                            p = sn.think_params{h};
                            P{h}.setThinkTime(jline.lang.processes.Pareto(p(1), p(2)));
                        case ProcessType.WEIBULL
                            p = sn.think_params{h};
                            P{h}.setThinkTime(jline.lang.processes.Weibull(p(1), p(2)));
                        case ProcessType.LOGNORMAL
                            p = sn.think_params{h};
                            P{h}.setThinkTime(jline.lang.processes.Lognormal(p(1), p(2)));
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

            %% Loop precedences (POST_LOOP)
            % sn.graph loop encoding -- see _kb/04-networkstruct.md
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
                jdist = javaObject('jline.lang.processes.Exp', line_dist.getParam(1).paramValue);
            elseif isa(line_dist, "APH")
                alpha = line_dist.getParam(1).paramValue;
                T = line_dist.getParam(2).paramValue;
                jline_alpha = java.util.ArrayList();
                for i = 1:length(alpha)
                    jline_alpha.add(alpha(i));
                end
                jline_T = JLINE.from_line_matrix(T);
                jdist = javaObject('jline.lang.processes.APH', jline_alpha, jline_T);
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
                jdist = javaObject('jline.lang.processes.Coxian', jline_mu, jline_phi);
            elseif isa(line_dist, 'Det')
                jdist = javaObject('jline.lang.processes.Det', line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'DiscreteSampler')
                popularity_p = JLINE.from_line_matrix(line_dist.getParam(1).paramValue);
                popularity_val = JLINE.from_line_matrix(line_dist.getParam(2).paramValue);
                jdist = javaObject('jline.lang.processes.DiscreteSampler', popularity_p, popularity_val);
            elseif isa(line_dist, 'Erlang')
                jdist = javaObject('jline.lang.processes.Erlang', line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Gamma')
                jdist = javaObject('jline.lang.processes.Gamma', line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, "HyperExp")
                jdist = javaObject('jline.lang.processes.HyperExp', line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue, line_dist.getParam(3).paramValue);
            elseif isa(line_dist, 'Lognormal')
                jdist = javaObject('jline.lang.processes.Lognormal', line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Pareto')
                jdist = javaObject('jline.lang.processes.Pareto', line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'MAP')
                D0 = line_dist.D(0);
                D1 = line_dist.D(1);
                jdist = javaObject('jline.lang.processes.MAP', JLINE.from_line_matrix(D0), JLINE.from_line_matrix(D1));
            elseif isa(line_dist, 'MMPP2')
                lambda0 =  line_dist.getParam(1).paramValue;
                lambda1 =  line_dist.getParam(2).paramValue;
                sigma0 =  line_dist.getParam(3).paramValue;
                sigma1 =  line_dist.getParam(4).paramValue;
                jdist = javaObject('jline.lang.processes.MMPP2', lambda0, lambda1, sigma0, sigma1);
            elseif isa(line_dist, 'NHPP')
                jdist = javaObject('jline.lang.processes.NHPP', line_dist.getBreakpoints(), line_dist.getRates(), logical(line_dist.isCyclic()));
            elseif isa(line_dist, 'MAPt') || isa(line_dist, 'PHt')
                % Both sides hold the schedule as per-segment matrix lists; the
                % JAR constructors take java.util.List<Matrix>
                if isa(line_dist, 'MAPt')
                    segA = line_dist.getD0Segments();
                    segB = line_dist.getD1Segments();
                else
                    segA = line_dist.getAlphaSegments();
                    segB = line_dist.getSSegments();
                end
                jA = javaObject('java.util.ArrayList');
                jB = javaObject('java.util.ArrayList');
                for k = 1:numel(segA)
                    jA.add(JLINE.from_line_matrix(segA{k}));
                    jB.add(JLINE.from_line_matrix(segB{k}));
                end
                jdist = javaObject(['jline.lang.processes.' class(line_dist)], ...
                    line_dist.getBreakpoints(), jA, jB, logical(line_dist.isCyclic()));
            elseif isa(line_dist, 'BMAP') % before MarkedMAP (BMAP < MarkedMAP)
                % Both sides use MarkedMAP layout {D0,D1_total,D1..DK} -- see _kb/12-interfaces-and-docs.md
                nmp = length(line_dist.process);
                jD = javaArray('jline.util.matrix.Matrix', nmp);
                for k = 1:nmp
                    jD(k) = JLINE.from_line_matrix(line_dist.process{k});
                end
                jdist = javaObject('jline.lang.processes.BMAP', javaObject('jline.util.matrix.MatrixCell', jD));
            elseif isa(line_dist, 'MarkedMMPP')
                nmp = length(line_dist.params); % {D0, D1, D11..D1K}
                jD = javaArray('jline.util.matrix.Matrix', nmp);
                for k = 1:nmp
                    jD(k) = JLINE.from_line_matrix(line_dist.getParam(k).paramValue);
                end
                jdist = javaObject('jline.lang.processes.MarkedMMPP', javaObject('jline.util.matrix.MatrixCell', jD));
            elseif isa(line_dist, 'MarkedMAP')
                nmp = length(line_dist.params); % {D0, D1, D11..D1K}
                jD = javaArray('jline.util.matrix.Matrix', nmp);
                for k = 1:nmp
                    jD(k) = JLINE.from_line_matrix(line_dist.getParam(k).paramValue);
                end
                jdist = javaObject('jline.lang.processes.MarkedMAP', javaObject('jline.util.matrix.MatrixCell', jD));
            elseif isa(line_dist, 'DMAP')
                jdist = javaObject('jline.lang.processes.DMAP', javaObject('jline.util.matrix.MatrixCell', JLINE.from_line_matrix(line_dist.getParam(1).paramValue), JLINE.from_line_matrix(line_dist.getParam(2).paramValue)));
            elseif isa(line_dist, 'ME')
                jdist = javaObject('jline.lang.processes.ME', JLINE.from_line_matrix(line_dist.getParam(1).paramValue), JLINE.from_line_matrix(line_dist.getParam(2).paramValue));
            elseif isa(line_dist, 'RAP')
                jdist = javaObject('jline.lang.processes.RAP', JLINE.from_line_matrix(line_dist.getParam(1).paramValue), JLINE.from_line_matrix(line_dist.getParam(2).paramValue));
            elseif isa(line_dist, 'MMDP2') % before MMDP (MMDP2 < MMDP)
                jdist = javaObject('jline.lang.processes.MMDP2', line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue, line_dist.getParam(3).paramValue, line_dist.getParam(4).paramValue);
            elseif isa(line_dist, 'MMDP')
                jdist = javaObject('jline.lang.processes.MMDP', JLINE.from_line_matrix(line_dist.getParam(1).paramValue), JLINE.from_line_matrix(line_dist.getParam(2).paramValue));
            elseif isa(line_dist, 'Normal')
                jdist = javaObject('jline.lang.processes.Normal', line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Bernoulli')
                jdist = javaObject('jline.lang.processes.Bernoulli', line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'Geometric')
                jdist = javaObject('jline.lang.processes.Geometric', line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'Poisson')
                jdist = javaObject('jline.lang.processes.Poisson', line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'Binomial')
                jdist = javaObject('jline.lang.processes.Binomial', int32(line_dist.getParam(1).paramValue), line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'DiscreteUniform')
                jdist = javaObject('jline.lang.processes.DiscreteUniform', line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'EmpiricalCDF')
                % MATLAB stores data as [cdf, x] column pairs; the JAR two-argument
                % constructor takes (cdfdata, xdata)
                if size(line_dist.data, 2) >= 2
                    jdist = javaObject('jline.lang.processes.EmpiricalCDF', JLINE.from_line_matrix(line_dist.data(:,1)), JLINE.from_line_matrix(line_dist.data(:,2)));
                else
                    jdist = javaObject('jline.lang.processes.EmpiricalCDF', JLINE.from_line_matrix(line_dist.data));
                end
            elseif isa(line_dist, 'PH')
                alpha = line_dist.getParam(1).paramValue;
                T = line_dist.getParam(2).paramValue;
                jdist = javaObject('jline.lang.processes.PH', JLINE.from_line_matrix(alpha), JLINE.from_line_matrix(T));
            elseif isa(line_dist, 'Uniform')
                jdist = javaObject('jline.lang.processes.Uniform', line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Weibull')
                jdist = javaObject('jline.lang.processes.Weibull', line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Zipf')
                jdist = javaObject('jline.lang.processes.Zipf', line_dist.getParam(3).paramValue, line_dist.getParam(4).paramValue);
            elseif isa(line_dist, 'Immediate')
                jdist = javaObject('jline.lang.processes.Immediate');
            elseif isempty(line_dist) || isa(line_dist, 'Disabled')
                jdist = javaObject('jline.lang.processes.Disabled');
                return;
            elseif isa(line_dist, 'Trace') % before Replayer (Trace < Replayer)
                jdist = javaObject('jline.lang.processes.Trace', line_dist.params{1}.paramValue);
            elseif isa(line_dist, 'Replayer')
                jdist = javaObject('jline.lang.processes.Replayer', line_dist.params{1}.paramValue);
            elseif isa(line_dist, 'Prior')
                % Convert Prior: each alternative is a distribution
                dists = line_dist.getParam(1).paramValue;
                probs = line_dist.getParam(2).paramValue;
                jdists = java.util.ArrayList();
                for k = 1:length(dists)
                    jdists.add(JLINE.from_line_distribution(dists{k}));
                end
                jdist = javaObject('jline.lang.processes.Prior', jdists, probs);
            else
                line_error(mfilename,'Distribution not supported by JLINE.');
            end
        end

        function [jRemDist, jRemPol] = from_line_signal_removal(line_class)
            % FROM_LINE_SIGNAL_REMOVAL Marshal a signal class's batch-removal
            % distribution and removal policy. Returns empties when the class
            % uses the defaults (remove exactly 1, RANDOM policy), so callers
            % can keep using the short JAR constructors in that case.
            jRemDist = [];
            jRemPol = [];
            if isprop(line_class, 'removalDistribution') && ~isempty(line_class.removalDistribution) ...
                    && ~isa(line_class.removalDistribution, 'Disabled')
                jRemDist = JLINE.from_line_distribution(line_class.removalDistribution);
            end
            if isprop(line_class, 'removalPolicy') && ~isempty(line_class.removalPolicy)
                if ~isempty(jRemDist) || line_class.removalPolicy ~= RemovalPolicy.RANDOM
                    jRemPol = jline.lang.constant.RemovalPolicy.fromID(int32(line_class.removalPolicy));
                end
            elseif ~isempty(jRemDist)
                jRemPol = jline.lang.constant.RemovalPolicy.RANDOM;
            end
        end

        function matlab_dist = from_jline_distribution(jdist)
            if isa(jdist, 'jline.lang.processes.Exp')
                matlab_dist = Exp(jdist.getRate());
            elseif isa(jdist, 'jline.lang.processes.Det')
                matlab_dist = Det(jdist.getParam(1).getValue);
            elseif isa(jdist, 'jline.lang.processes.Erlang')
                matlab_dist = Erlang(jdist.getParam(1).getValue(),jdist.getNumberOfPhases());
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
            elseif isa(jdist, 'jline.lang.processes.MAP')
                D0 = JLINE.from_jline_matrix(jdist.D(0));
                D1 = JLINE.from_jline_matrix(jdist.D(1));
                matlab_dist = MAP({D0, D1});
            elseif isa(jdist, 'jline.lang.processes.APH')
                alpha = JLINE.from_jline_matrix(jdist.getInitProb());
                T = JLINE.from_jline_matrix(jdist.getSubgenerator());
                matlab_dist = APH(alpha(:)', T);
            elseif isa(jdist, 'jline.lang.processes.PH')
                alpha = JLINE.from_jline_matrix(jdist.getInitProb());
                T = JLINE.from_jline_matrix(jdist.getSubgenerator());
                matlab_dist = PH(alpha(:)', T);
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
            elseif isa(jdist, 'jline.lang.processes.Zipf')
                matlab_dist = Zipf(jdist.getParam(3).getValue, jdist.getParam(4).getValue);
            elseif isa(jdist, 'jline.lang.processes.DiscreteSampler')
                jpMat = jdist.getParam(1).getValue;
                jxMat = jdist.getParam(2).getValue;
                p = zeros(1, jpMat.length);
                x = zeros(1, jxMat.length);
                for i = 1:jpMat.length
                    p(i) = jpMat.get(i-1);
                end
                for i = 1:jxMat.length
                    x(i) = jxMat.get(i-1);
                end
                matlab_dist = DiscreteSampler(p, x);
            elseif isa(jdist, 'jline.lang.processes.Immediate')
                matlab_dist = Immediate();
            elseif isa(jdist, 'jline.lang.processes.Disabled')
                matlab_dist = Disabled();
            elseif isa(jdist, 'jline.lang.processes.MMPP2')
                matlab_dist = MMPP2(jdist.getParam(1).getValue, jdist.getParam(2).getValue, jdist.getParam(3).getValue, jdist.getParam(4).getValue);
            elseif isa(jdist, 'jline.lang.processes.NHPP')
                jBp = jdist.getBreakpoints();
                jRates = jdist.getRates();
                bp = zeros(1, length(jBp));
                rates = zeros(1, length(jRates));
                for k = 1:length(jBp)
                    bp(k) = jBp(k);
                end
                for k = 1:length(jRates)
                    rates(k) = jRates(k);
                end
                matlab_dist = NHPP(bp, rates, logical(jdist.isCyclic()));
            elseif isa(jdist, 'jline.lang.processes.MAPt') || isa(jdist, 'jline.lang.processes.PHt')
                jBp = jdist.getBreakpoints();
                bp = zeros(1, length(jBp));
                for k = 1:length(jBp)
                    bp(k) = jBp(k);
                end
                isMAPt = isa(jdist, 'jline.lang.processes.MAPt');
                if isMAPt
                    jA = jdist.getD0Segments(); jB = jdist.getD1Segments();
                else
                    jA = jdist.getAlphaSegments(); jB = jdist.getSSegments();
                end
                n = jA.size();
                segA = cell(1, n); segB = cell(1, n);
                for k = 1:n
                    segA{k} = JLINE.from_jline_matrix(jA.get(k-1));
                    segB{k} = JLINE.from_jline_matrix(jB.get(k-1));
                end
                if isMAPt
                    matlab_dist = MAPt(bp, segA, segB, logical(jdist.isCyclic()));
                else
                    matlab_dist = PHt(bp, segA, segB, logical(jdist.isCyclic()));
                end
            elseif isa(jdist, 'jline.lang.processes.Prior')
                % Convert Prior from JAR to MATLAB
                jdists = jdist.getDistributions();
                nalt = jdists.size();
                dists = cell(1, nalt);
                for k = 1:nalt
                    dists{k} = JLINE.from_jline_distribution(jdists.get(k-1));
                end
                probs = jdist.getProbabilities();
                matlab_dist = Prior(dists, probs);
            elseif isa(jdist, 'jline.lang.processes.Normal')
                matlab_dist = Normal(jdist.getParam(1).getValue, jdist.getParam(2).getValue);
            elseif isa(jdist, 'jline.lang.processes.Bernoulli')
                matlab_dist = Bernoulli(jdist.getParam(1).getValue);
            elseif isa(jdist, 'jline.lang.processes.Geometric')
                matlab_dist = Geometric(jdist.getParam(1).getValue);
            elseif isa(jdist, 'jline.lang.processes.Poisson')
                matlab_dist = Poisson(jdist.getParam(1).getValue);
            elseif isa(jdist, 'jline.lang.processes.Binomial')
                matlab_dist = Binomial(double(jdist.getParam(1).getValue), jdist.getParam(2).getValue);
            elseif isa(jdist, 'jline.lang.processes.DiscreteUniform')
                matlab_dist = DiscreteUniform(jdist.getParam(1).getValue, jdist.getParam(2).getValue);
            elseif isa(jdist, 'jline.lang.processes.EmpiricalCDF')
                d = JLINE.from_jline_matrix(jdist.getData());
                if size(d, 2) >= 2 % stored as [cdf, x]; MATLAB ctor is (xdata, cdfdata)
                    matlab_dist = EmpiricalCDF(d(:,2), d(:,1));
                else
                    matlab_dist = EmpiricalCDF(d);
                end
            elseif isa(jdist, 'jline.lang.processes.Trace') % before Replayer (Trace < Replayer)
                v = jdist.getParam(1).getValue;
                if isa(v, 'java.lang.String'), v = char(v); end
                matlab_dist = Trace(v);
            elseif isa(jdist, 'jline.lang.processes.Replayer')
                v = jdist.getParam(1).getValue;
                if isa(v, 'java.lang.String'), v = char(v); end
                matlab_dist = Replayer(v);
            elseif isa(jdist, 'jline.lang.processes.DMAP')
                matlab_dist = DMAP(JLINE.from_jline_matrix(jdist.getParam(1).getValue), JLINE.from_jline_matrix(jdist.getParam(2).getValue));
            elseif isa(jdist, 'jline.lang.processes.ME')
                matlab_dist = ME(JLINE.from_jline_matrix(jdist.getParam(1).getValue), JLINE.from_jline_matrix(jdist.getParam(2).getValue));
            elseif isa(jdist, 'jline.lang.processes.RAP')
                matlab_dist = RAP(JLINE.from_jline_matrix(jdist.getParam(1).getValue), JLINE.from_jline_matrix(jdist.getParam(2).getValue));
            elseif isa(jdist, 'jline.lang.processes.MMDP2') % before MMDP (MMDP2 < MMDP)
                matlab_dist = MMDP2(jdist.getParam(1).getValue, jdist.getParam(2).getValue, jdist.getParam(3).getValue, jdist.getParam(4).getValue);
            elseif isa(jdist, 'jline.lang.processes.MMDP')
                matlab_dist = MMDP(JLINE.from_jline_matrix(jdist.getParam(1).getValue), JLINE.from_jline_matrix(jdist.getParam(2).getValue));
            elseif isa(jdist, 'jline.lang.processes.BMAP') % before MarkedMAP (BMAP < MarkedMAP)
                % JAR cell is MarkedMAP layout; rebuilt to standard {D0,D1,...} -- see _kb/12-interfaces-and-docs.md
                proc = jdist.getProcess();
                K = proc.size() - 2;
                D = cell(1, K+1);
                D{1} = JLINE.from_jline_matrix(proc.get(0));
                for k = 1:K
                    D{1+k} = JLINE.from_jline_matrix(proc.get(1+k));
                end
                matlab_dist = BMAP(D);
            elseif isa(jdist, 'jline.lang.processes.MarkedMMPP')
                proc = jdist.getProcess(); % {D0, D1, D11..D1K}
                K = proc.size() - 2;
                D = cell(1, K+2);
                for k = 1:K+2
                    D{k} = JLINE.from_jline_matrix(proc.get(k-1));
                end
                matlab_dist = MarkedMMPP(D, K);
            elseif isa(jdist, 'jline.lang.processes.MarkedMAP')
                proc = jdist.getProcess(); % {D0, D1, D11..D1K}
                K = proc.size() - 2;
                D = cell(1, K+2);
                for k = 1:K+2
                    D{k} = JLINE.from_jline_matrix(proc.get(k-1));
                end
                matlab_dist = MarkedMAP(D, K);
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
            if (isa(line_node,'Sink')) || isa(line_node, 'ClassSwitch') || isa(line_node, 'Fork') || isa(line_node, 'Join') || isa(line_node, 'Place') || isa(line_node, 'Transition') || isa(line_node, 'Cache')
                return;
            end
            for n = 1:job_classes.size()
                if (isa(line_node, 'Queue') || isa(line_node, 'Delay'))
                    jdist = jline_node.getServiceProcess(job_classes.get(n-1));
                    matlab_dist = JLINE.from_jline_distribution(jdist);
                    weight = jline_node.getSchedStrategyPar(job_classes.get(n-1));
                    line_node.setService(line_classes{n}, matlab_dist, weight);
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
                node_object = javaObject('jline.lang.nodes.Delay', jnetwork, line_node.getName);
            elseif isa(line_node, 'Queue')
                switch line_node.schedStrategy
                    case SchedStrategy.INF
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.INF);
                    case SchedStrategy.FCFS
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFS);
                    case SchedStrategy.LCFS
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFS);
                    case SchedStrategy.SIRO
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.SIRO);
                    case SchedStrategy.SJF
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.SJF);
                    case SchedStrategy.LJF
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LJF);
                    case SchedStrategy.PS
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.PS);
                    case SchedStrategy.DPS
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.DPS);
                    case SchedStrategy.GPS
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.GPS);
                    case SchedStrategy.SEPT
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.SEPT);
                    case SchedStrategy.LEPT
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LEPT);
                    case SchedStrategy.HOL
                        % HOL maps to the JAR's own HOL, not FCFSPRIO -- see _kb/12-interfaces-and-docs.md
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.HOL);
                    case SchedStrategy.FCFSPRIO
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFSPRIO);
                    case SchedStrategy.FORK
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FORK);
                    case SchedStrategy.EXT
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.EXT);
                    case SchedStrategy.REF
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.REF);
                    case SchedStrategy.LCFSPR
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFSPR);
                    case SchedStrategy.LCFSPI
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFSPI);
                    case SchedStrategy.LCFSPRIO
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFSPRIO);
                    case SchedStrategy.LCFSPRPRIO
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFSPRPRIO);
                    case SchedStrategy.LCFSPIPRIO
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LCFSPIPRIO);
                    case SchedStrategy.FCFSPR
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFSPR);
                    case SchedStrategy.FCFSPI
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFSPI);
                    case SchedStrategy.FCFSPRPRIO
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFSPRPRIO);
                    case SchedStrategy.FCFSPIPRIO
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FCFSPIPRIO);
                    case SchedStrategy.PSPRIO
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.PSPRIO);
                    case SchedStrategy.DPSPRIO
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.DPSPRIO);
                    case SchedStrategy.GPSPRIO
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.GPSPRIO);
                    case SchedStrategy.POLLING
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.POLLING);
                    case SchedStrategy.SRPT
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.SRPT);
                    case SchedStrategy.SRPTPRIO
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.SRPTPRIO);
                    case SchedStrategy.PSJF
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.PSJF);
                    case SchedStrategy.FB
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FB);
                    case SchedStrategy.LRPT
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LRPT);
                    case SchedStrategy.EDD
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.EDD);
                    case SchedStrategy.EDF
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.EDF);
                    case SchedStrategy.LPS
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.LPS);
                    case SchedStrategy.SETF
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.SETF);
                    case SchedStrategy.FSP
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.FSP);
                    case SchedStrategy.PAS
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.PAS);
                    case SchedStrategy.OI
                        node_object = javaObject('jline.lang.nodes.Queue', jnetwork, line_node.getName, jline.lang.constant.SchedStrategy.OI);
                    otherwise
                        line_error(mfilename, sprintf('JLINE conversion does not support the %s scheduling strategy yet.', char(SchedStrategy.toText(line_node.schedStrategy))));
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
                    cdPeakJava = line_node.lcdScalingPeak;
                    if isscalar(cdPeakJava)
                        cdPeakJava = repmat(cdPeakJava, 1, sn.nclasses);
                    end
                    node_object.setLimitedClassDependence(JLINE.handle_to_serializablefun(line_node.lcdScaling, sn), JLINE.from_line_matrix(cdPeakJava(:)'));
                end
                if ~isempty(line_node.ljdScaling)
                    if isempty(sn)
                        line_error(mfilename, "Joint-dependent models require sn struct for MATLAB-to-JAVA translation.");
                    end
                    jdPeakJava = line_node.ljdScalingPeak;
                    if isscalar(jdPeakJava)
                        jdPeakJava = repmat(jdPeakJava, 1, sn.nclasses);
                    end
                    node_object.setLimitedJointDependence(JLINE.handle_to_serializablefun(line_node.ljdScaling, sn), JLINE.from_line_matrix(jdPeakJava(:)'));
                end
                % Set queue capacity if finite
                if ~isinf(line_node.cap)
                    node_object.setCapacity(line_node.cap);
                end
                % Transfer LPS job limit (stored in schedStrategyPar(1))
                if line_node.schedStrategy == SchedStrategy.LPS && ~isempty(line_node.schedStrategyPar) && line_node.schedStrategyPar(1) >= 1
                    node_object.setLimit(int32(line_node.schedStrategyPar(1)));
                end
            elseif isa(line_node, 'Source')
                node_object = javaObject('jline.lang.nodes.Source', jnetwork, line_node.getName);
            elseif isa(line_node, 'Sink')
                node_object = javaObject('jline.lang.nodes.Sink', jnetwork, line_node.getName);
            elseif isa(line_node, 'Router')
                node_object = javaObject('jline.lang.nodes.Router', jnetwork, line_node.getName);
            elseif isa(line_node, 'ClassSwitch')
                node_object = javaObject('jline.lang.nodes.ClassSwitch', jnetwork, line_node.getName);
            elseif isa(line_node, 'Fork')
                node_object = javaObject('jline.lang.nodes.Fork', jnetwork, line_node.name);
                node_object.setTasksPerLink(line_node.output.tasksPerLink);
            elseif isa(line_node, 'Join')
                node_object = javaObject('jline.lang.nodes.Join', jnetwork, line_node.name, forkNode);
            elseif isa(line_node, 'Logger')
                node_object = javaObject('jline.lang.nodes.Logger', jnetwork, line_node.name, [line_node.filePath,line_node.fileName]);
                % Output-field flags transferred unconditionally (JAR defaults all false) -- see _kb/12-interfaces-and-docs.md
                node_object.setStartTime(strcmpi(char(line_node.getStartTime), 'true'));
                node_object.setLoggerName(strcmpi(char(line_node.getLoggerName), 'true'));
                node_object.setTimestamp(strcmpi(char(line_node.getTimestamp), 'true'));
                node_object.setJobID(strcmpi(char(line_node.getJobID), 'true'));
                node_object.setJobClass(strcmpi(char(line_node.getJobClass), 'true'));
                node_object.setTimeSameClass(strcmpi(char(line_node.getTimeSameClass), 'true'));
                node_object.setTimeAnyClass(strcmpi(char(line_node.getTimeAnyClass), 'true'));
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
                    % graph is a per-item cell array of (h+1)x(h+1) matrices
                    gcells = line_node.graph;
                    if ~iscell(gcells), gcells = {gcells}; end
                    jGraph = javaArray('jline.util.matrix.Matrix', numel(gcells));
                    for gidx = 1:numel(gcells)
                        jGraph(gidx) = JLINE.from_line_matrix(gcells{gidx});
                    end
                    node_object = javaObject('jline.lang.nodes.Cache', jnetwork, line_node.name, nitems, JLINE.from_line_matrix(line_node.itemLevelCap), repStrategy, jGraph);
                else
                    node_object = javaObject('jline.lang.nodes.Cache', jnetwork, line_node.name, nitems, JLINE.from_line_matrix(line_node.itemLevelCap), repStrategy);
                end
                if isprop(line_node,'admissionProb') && ~isempty(line_node.admissionProb)
                    node_object.setAdmissionProb(line_node.admissionProb);
                end
                % per-item storage costs and per-list cost caps (ton21cache Sec. IX)
                if isprop(line_node,'itemSize') && ~isempty(line_node.itemSize)
                    node_object.setItemSizes(JLINE.from_line_matrix(line_node.itemSize(:)'));
                end
                if isprop(line_node,'costCap') && ~isempty(line_node.costCap)
                    node_object.setCostCaps(JLINE.from_line_matrix(line_node.costCap(:)'));
                end
            elseif isa(line_node, 'Place')
                if line_node.isQueueing()
                    % Reconstructed with its scheduling strategy for setService -- see _kb/12-interfaces-and-docs.md
                    jsched = jline.lang.constant.SchedStrategy.fromText(SchedStrategy.toText(line_node.schedStrategy));
                    node_object = javaObject('jline.lang.nodes.Place', jnetwork, line_node.getName, jsched);
                else
                    node_object = javaObject('jline.lang.nodes.Place', jnetwork, line_node.getName);
                end
            elseif isa(line_node, 'Transition')
                % A MATLAB function handle cannot cross into a Java
                % SerializableFunction, so g(marking) would be dropped and the
                % solve would return nominal (unscaled) firing rates silently
                if isprop(line_node, 'firingRateDependence') && ...
                        any(~cellfun(@isempty, line_node.firingRateDependence))
                    line_error(mfilename, sprintf(['Transition ''%s'' has a marking-dependent firing rate, ' ...
                        'which cannot be marshalled to JLINE. Use lang=''matlab'' or lang=''python''.'], ...
                        line_node.getName));
                end
                node_object = javaObject('jline.lang.nodes.Transition', jnetwork, line_node.getName);
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
                    case 'PSPRIO'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.PSPRIO);
                    case 'DPSPRIO'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.DPSPRIO);
                    case 'GPSPRIO'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.GPSPRIO);
                    case 'LCFSPI'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.LCFSPI);
                    case 'LCFSPRIO'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.LCFSPRIO);
                    case 'LCFSPRPRIO'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.LCFSPRPRIO);
                    case 'LCFSPIPRIO'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.LCFSPIPRIO);
                    case 'FCFSPR'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.FCFSPR);
                    case 'FCFSPI'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.FCFSPI);
                    case 'FCFSPRPRIO'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.FCFSPRPRIO);
                    case 'FCFSPIPRIO'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.FCFSPIPRIO);
                    case 'POLLING'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.POLLING);
                    case 'EDD'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.EDD);
                    case 'EDF'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.EDF);
                    case 'LPS'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.LPS);
                    case 'SETF'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.SETF);
                    case 'FSP'
                        node_object = Queue(model, jline_node.getName.toCharArray', SchedStrategy.FSP);
                    case 'PAS'
                        % PAS rate function is a Java SerializableFunction, not recoverable as a MATLAB handle
                        line_error(mfilename, 'JLINE-to-LINE conversion does not support PAS queues (service rate function not recoverable).');
                    otherwise
                        line_error(mfilename, sprintf('JLINE-to-LINE conversion does not support the %s scheduling strategy yet.', char(schedStrategy.name())));
                end
                node_object.setNumberOfServers(jline_node.getNumberOfServers);
                cap = jline_node.getCap();
                if cap < intmax && cap > 0
                    node_object.setCapacity(cap);
                end
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
            elseif isa(jline_node, 'jline.lang.nodes.Cache')
                numItems = jline_node.getNumberOfItems();
                itemLevelCap = JLINE.from_jline_matrix(jline_node.getItemLevelCap());
                replPolicy = jline_node.getReplacementStrategy();
                switch char(replPolicy)
                    case 'LRU'
                        rp = ReplacementStrategy.LRU;
                    case 'FIFO'
                        rp = ReplacementStrategy.FIFO;
                    case 'RR'
                        rp = ReplacementStrategy.RR;
                    otherwise
                        rp = ReplacementStrategy.LRU;
                end
                node_object = Cache(model, jline_node.getName.toCharArray', numItems, itemLevelCap, rp);
                % hitClass, missClass, and popularity set later in jline_to_line
            elseif isa(jline_node, 'jline.lang.nodes.Fork')
                node_object = Fork(model, jline_node.getName.toCharArray');
                tpl = jline_node.getOutput().tasksPerLink;
                if tpl > 1
                    node_object.setTasksPerLink(tpl);
                end
            elseif isa(jline_node, 'jline.lang.nodes.Join')
                node_object = Join(model, jline_node.getName.toCharArray');
                % joinOf is set later in jline_to_line after all nodes are created
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
                [jRemDist, jRemPol] = JLINE.from_line_signal_removal(line_class);
                if isempty(jRemDist) && isempty(jRemPol)
                    node_class = javaObject('jline.lang.ClosedSignal', jnetwork, line_class.getName, jSignalType, jnetwork.getNodeByName(line_class.refstat.getName), line_class.priority);
                else
                    node_class = javaObject('jline.lang.ClosedSignal', jnetwork, line_class.getName, jSignalType, jnetwork.getNodeByName(line_class.refstat.getName), line_class.priority, jRemDist, jRemPol);
                end
            elseif isa(line_class, 'Signal') || isa(line_class, 'OpenSignal')
                % Signal/OpenSignal -> jline.lang.Signal (includes CATASTROPHE type)
                jSignalType = jline.lang.constant.SignalType.fromID(line_class.signalType);
                [jRemDist, jRemPol] = JLINE.from_line_signal_removal(line_class);
                if isempty(jRemDist) && isempty(jRemPol)
                    node_class = javaObject('jline.lang.Signal', jnetwork, line_class.getName, jSignalType, line_class.priority);
                else
                    node_class = javaObject('jline.lang.Signal', jnetwork, line_class.getName, jSignalType, line_class.priority, jRemDist, jRemPol);
                end
            elseif isa(line_class, 'OpenClass')
                node_class = javaObject('jline.lang.OpenClass', jnetwork, line_class.getName, line_class.priority);
            elseif isa(line_class, 'SelfLoopingClass')
                node_class = javaObject('jline.lang.SelfLoopingClass', jnetwork, line_class.getName, line_class.population, jnetwork.getNodeByName(line_class.refstat.getName), line_class.priority);
            elseif isa(line_class, 'ClosedClass')
                node_class = javaObject('jline.lang.ClosedClass', jnetwork, line_class.getName, line_class.population, jnetwork.getNodeByName(line_class.refstat.getName), line_class.priority);
            else
                line_error(mfilename,'Class type not supported by JLINE.');
            end
            % Transfer relative deadline (used by EDD/EDF scheduling)
            if isprop(line_class, 'deadline') && ~isempty(line_class.deadline) && ~isinf(line_class.deadline)
                node_class.setDeadline(line_class.deadline);
            end
            if line_class.isReferenceClass()
                node_class.setReferenceClass(true);
            end
        end

        function node_class = from_jline_class(jclass, model)
            % Check signal classes first (their base classes would match below)
            if isa(jclass, 'jline.lang.ClosedSignal')
                [remDist, remPol] = JLINE.from_jline_signal_removal(jclass);
                node_class = ClosedSignal(model, jclass.getName.toCharArray', jclass.getSignalType().getID(), model.getNodeByName(jclass.getReferenceStation.getName), jclass.getPriority, remDist, remPol);
            elseif isa(jclass, 'jline.lang.OpenSignal')
                node_class = OpenSignal(model, jclass.getName.toCharArray', jclass.getSignalType().getID(), jclass.getPriority);
            elseif isa(jclass, 'jline.lang.Signal')
                [remDist, remPol] = JLINE.from_jline_signal_removal(jclass);
                node_class = Signal(model, jclass.getName.toCharArray', jclass.getSignalType().getID(), jclass.getPriority, remDist, remPol);
            elseif isa(jclass, 'jline.lang.OpenClass')
                node_class = OpenClass(model, jclass.getName.toCharArray', jclass.getPriority);
            elseif isa(jclass, 'jline.lang.SelfLoopingClass')
                node_class = SelfLoopingClass(model, jclass.getName.toCharArray', jclass.getNumberOfJobs, model.getNodeByName(jclass.getReferenceStation.getName), jclass.getPriority);
            elseif isa(jclass, 'jline.lang.ClosedClass')
                node_class = ClosedClass(model, jclass.getName.toCharArray', jclass.getNumberOfJobs, model.getNodeByName(jclass.getReferenceStation.getName), jclass.getPriority);
            else
                line_error(mfilename,'Class type not supported by JLINE.');
            end
            % Transfer relative deadline (used by EDD/EDF scheduling)
            dl = jclass.getDeadline();
            if isfinite(dl)
                node_class.deadline = dl;
            end
        end

        function [remDist, remPol] = from_jline_signal_removal(jclass)
            % FROM_JLINE_SIGNAL_REMOVAL Read back a JAR signal's batch-removal
            % distribution and policy for the MATLAB signal constructors.
            jRemDist = jclass.getRemovalDistribution();
            if isempty(jRemDist)
                remDist = [];
            else
                remDist = JLINE.from_jline_distribution(jRemDist);
            end
            jRemPol = jclass.getRemovalPolicy();
            if isempty(jRemPol)
                remPol = RemovalPolicy.RANDOM;
            else
                remPol = jRemPol.getID();
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
                % sn.rtorig already excludes auto-added ClassSwitch nodes (station-indexed) -- see _kb/12-interfaces-and-docs.md
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
                                if useLinkMethod
                                    % No routing-matrix entries for RAND under useLinkMethod -- see _kb/12-interfaces-and-docs.md
                                else
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
                            case RoutingStrategy.SDR
                                % Krzesinski state-dependent routing: the JAR
                                % takes the whole subnetwork in one call, so the
                                % declared branch topology is translated from
                                % MATLAB node objects to the JAR's own nodes.
                                % See _kb/16-state-dependent-routing.md
                                outlinks_i=find(connections(i,:))';
                                if useLinkMethod
                                    line_error(mfilename,'State-dependent routing cannot be used together with the link() command.');
                                end
                                for j= outlinks_i(:)'
                                    jdest_idx = matlab2java_node_idx(j);
                                    if jdest_idx >= 0
                                        jmodel.addLink(jnodes.get(jnode_idx), jnodes.get(jdest_idx));
                                    end
                                end
                                decl = output_strat{3}{1};
                                B = numel(decl.branch);
                                jbranches = java.util.ArrayList();
                                jbranches.add(java.util.ArrayList());  % index 1 is the complement M-V
                                for b = 2:B
                                    jb = java.util.ArrayList();
                                    for q = 1:numel(decl.branch{b})
                                        jb.add(jmodel.getNodeByName(decl.branch{b}{q}.getName()));
                                    end
                                    jbranches.add(jb);
                                end
                                jlevel = int32(decl.level(:)');
                                jC = decl.C(:)';
                                jd = decl.d;
                                jnodes.get(jnode_idx).setStateDepRouting(jclasses.get(k-1), ...
                                    jmodel.getNodeByName(decl.departure.getName()), ...
                                    jbranches, jlevel, jC, jd);
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
                            case RoutingStrategy.JSQ
                                jnodes.get(jnode_idx).setRouting(jclasses.get(k-1),jline.lang.constant.RoutingStrategy.JSQ);
                                outlinks_i=find(connections(i,:))';
                                if ~useLinkMethod
                                    for j= outlinks_i(:)'
                                        jdest_idx = matlab2java_node_idx(j);
                                        if jdest_idx >= 0
                                            jmodel.addLink(jnodes.get(jnode_idx), jnodes.get(jdest_idx));
                                        end
                                    end
                                end
                            case RoutingStrategy.SQ
                                jnodes.get(jnode_idx).setRouting(jclasses.get(k-1),jline.lang.constant.RoutingStrategy.SQ);
                                if length(output_strat) >= 3 && ~isempty(output_strat{3})
                                    dparam = output_strat{3}{1};
                                    jnodes.get(jnode_idx).setSQRouting(jclasses.get(k-1), int32(dparam));
                                end
                                outlinks_i=find(connections(i,:))';
                                if ~useLinkMethod
                                    for j= outlinks_i(:)'
                                        jdest_idx = matlab2java_node_idx(j);
                                        if jdest_idx >= 0
                                            jmodel.addLink(jnodes.get(jnode_idx), jnodes.get(jdest_idx));
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
            n_nodes = jnodes.size();
            network_nodes = model.getNodes;
            network_classes = model.getClasses;

            % Build name-to-MATLAB-node map (JAR and MATLAB may order nodes differently)
            node_by_name = configureDictionary('string','cell');
            for nn = 1:length(network_nodes)
                node_by_name{network_nodes{nn}.name} = network_nodes{nn};
            end

            connections = JLINE.from_jline_matrix(jnetwork.getConnectionMatrix());
            [row,col] = find(connections);
            for i=1:length(row)
                from_name = char(jnodes.get(row(i)-1).getName());
                to_name = char(jnodes.get(col(i)-1).getName());
                model.addLink(node_by_name{from_name}, node_by_name{to_name});
            end

            for n = 1 : n_nodes
                jnode = jnodes.get(n-1);
                cur_node = node_by_name{char(jnode.getName())};
                output_strategies = jnode.getOutputStrategies();
                n_strategies = output_strategies.size();
                for m = 1 : n_strategies
                    output_strat = output_strategies.get(m-1);
                    routing_strat = output_strat.getRoutingStrategy;
                    routing_strat_classidx = output_strat.getJobClass.getIndex();
                    switch char(routing_strat)
                        case 'RAND'
                            cur_node.setRouting(network_classes{routing_strat_classidx}, RoutingStrategy.RAND);
                        case 'RROBIN'
                            cur_node.setRouting(network_classes{routing_strat_classidx}, RoutingStrategy.RROBIN);
                        case 'WRROBIN'
                            dest = output_strat.getDestination();
                            if ~isempty(dest)
                                dest_name = char(dest.getName());
                                if isKey(node_by_name, dest_name)
                                    weight = output_strat.getProbability();
                                    cur_node.setRouting(network_classes{routing_strat_classidx}, RoutingStrategy.WRROBIN, node_by_name{dest_name}, weight);
                                end
                            end
                        case 'DISABLED'
                            cur_node.setRouting(network_classes{routing_strat_classidx}, RoutingStrategy.DISABLED);
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
            network_nodes = model.getNodes;

            % Build JAR-to-MATLAB node index mapping (node orders may differ)
            jar2ml = zeros(1, n_nodes);
            for jj = 1:n_nodes
                jar_name = char(jnodes.get(jj-1).getName());
                for mm = 1:length(network_nodes)
                    if strcmp(network_nodes{mm}.name, jar_name)
                        jar2ml(jj) = mm;
                        break;
                    end
                end
            end

            % rtorig is the matrix the JAR model was linked with, so it is
            % complete for EVERY class. The OutputStrategy scan below is not: a
            % class routed by a non-PROB strategy (RAND, RROBIN, ...) holds a
            % single entry with a NULL destination, so its whole row is invisible
            % there. Reading the strategies first therefore silently dropped
            % every connection used only by such a class -- a Delay self-loop
            % under RAND (cqn_mmpp2_service) and a Router's outgoing links
            % (cache_replc_routing), which severed the network. rtorig is
            % preferred whenever it is available and the scan is the fallback.
            % Its order is that of the nodes BEFORE link ran, i.e. without the
            % auto-added CS nodes; a mismatch means that broke -- see _kb/12.
            usedRtorig = false;
            sn = jnetwork.getStruct;
            n_ml_nodes = length(model.getNodes);
            if ~isempty(sn.rtorig)
                rtblocks = cell(n_classes, n_classes);
                sized = true;
                for r = 1:n_classes
                    for s = 1:n_classes
                        rtMat = JLINE.from_jline_matrix(sn.rtorig.get(jclasses.get(r-1)).get(jclasses.get(s-1)));
                        if ~isempty(rtMat)
                            if ~isequal(size(rtMat), [n_ml_nodes, n_ml_nodes])
                                sized = false;
                            end
                            rtblocks{r,s} = rtMat;
                            usedRtorig = usedRtorig || any(rtMat(:) > 0);
                        end
                    end
                end
                if sized && usedRtorig
                    for r = 1:n_classes
                        for s = 1:n_classes
                            if ~isempty(rtblocks{r,s})
                                P{r,s} = rtblocks{r,s};
                            end
                        end
                    end
                else
                    usedRtorig = false;
                end
            end

            if ~usedRtorig
                for n = 1 : n_nodes
                    jnode = jnodes.get(n-1);
                    output_strategies = jnode.getOutputStrategies();
                    n_strategies = output_strategies.size();
                    for m = 1 : n_strategies
                        output_strat = output_strategies.get(m-1);
                        dest = output_strat.getDestination();
                        if~isempty(dest) % disabled strategy
                            in_idx = jar2ml(jnetwork.getNodeIndex(jnode)+1);
                            out_idx = jar2ml(jnetwork.getNodeIndex(dest)+1);
                            if in_idx == 0 || out_idx == 0
                                continue; % node absent from the rebuilt model
                            end
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
            end

            model.link(P);

            % Restore non-PROB routing strategies (RROBIN, WRROBIN, etc.)
            % after link(), which sets all routing to PROB
            network_nodes = model.getNodes;
            network_classes = model.getClasses;
            node_by_name = configureDictionary('string','cell');
            for nn = 1:length(network_nodes)
                node_by_name{network_nodes{nn}.name} = network_nodes{nn};
            end
            for n = 1 : n_nodes
                jnode = jnodes.get(n-1);
                if ~isKey(node_by_name, char(jnode.getName()))
                    continue; % a node link() re-created carries link()'s routing
                end
                cur_node = node_by_name{char(jnode.getName())};
                output_strategies = jnode.getOutputStrategies();
                n_strategies = output_strategies.size();
                for m = 1 : n_strategies
                    output_strat = output_strategies.get(m-1);
                    routing_strat = output_strat.getRoutingStrategy;
                    routing_strat_classidx = output_strat.getJobClass.getIndex();
                    switch char(routing_strat)
                        % RAND and DISABLED are DECLARED, not derived: link(P)
                        % installs the resolved probabilities as PROB, so without
                        % this the round-tripped model reports a different
                        % strategy than the one it was built with (numerically
                        % equal only while the outlinks stay equiprobable).
                        case 'RAND'
                            cur_node.setRouting(network_classes{routing_strat_classidx}, RoutingStrategy.RAND);
                        case 'DISABLED'
                            cur_node.setRouting(network_classes{routing_strat_classidx}, RoutingStrategy.DISABLED);
                        case 'RROBIN'
                            cur_node.setRouting(network_classes{routing_strat_classidx}, RoutingStrategy.RROBIN);
                        case 'WRROBIN'
                            dest = output_strat.getDestination();
                            if ~isempty(dest)
                                dest_name = char(dest.getName());
                                if isKey(node_by_name, dest_name)
                                    % Clear stale PROB entries before first WRROBIN weight
                                    classIdx = routing_strat_classidx;
                                    if length(cur_node.output.outputStrategy) >= classIdx && ...
                                            length(cur_node.output.outputStrategy{1, classIdx}) >= 3
                                        curStrat = cur_node.output.outputStrategy{1, classIdx}{2};
                                        if ~strcmp(curStrat, 'WeightedRoundRobin')
                                            cur_node.output.outputStrategy{1, classIdx}{3} = {};
                                        end
                                    end
                                    weight = output_strat.getProbability();
                                    cur_node.setRouting(network_classes{routing_strat_classidx}, RoutingStrategy.WRROBIN, node_by_name{dest_name}, weight);
                                end
                            end
                    end
                end
            end

            % Invalidate cached struct after modifying routing strategies
            model.resetStruct();

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

            jnetwork = javaObject('jline.lang.Network', model.getName);
            % setChecks carried over before link() (SolverLN fast-mode layers) -- see _kb/12-interfaces-and-docs.md
            jnetwork.setChecks(logical(model.getChecks));
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
                                case DropStrategy.BBS
                                    jnodes{n}.setDropRule(jclasses{r}, jline.lang.constant.DropStrategy.BlockingBeforeService);
                                case DropStrategy.RSRD
                                    jnodes{n}.setDropRule(jclasses{r}, jline.lang.constant.DropStrategy.ReServiceOnRejection);
                                % WAITQ (-1, MATLAB's universal default) is never transferred -- see _kb/12-interfaces-and-docs.md
                            end
                        end
                    end
                    % Transfer per-class capacity limits (setChainCapacity/classCap)
                    if ~isempty(line_nodes{n}.classCap)
                        for r = 1:min(length(line_classes), length(line_nodes{n}.classCap))
                            if isfinite(line_nodes{n}.classCap(r)) && line_nodes{n}.classCap(r) >= 0
                                jnodes{n}.setClassCap(jclasses{r}, int32(line_nodes{n}.classCap(r)));
                            end
                        end
                    end
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
                            case PollingType.DECREMENTING
                                jPollingType = jline.lang.constant.PollingType.DECREMENTING;
                            otherwise
                                line_error(mfilename, sprintf('Unsupported polling type for the Java backend: %d.', PollingType.toId(pollingType)));
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

            % Transfer impatience features (reneging, balking, retrial) for queues
            % These require classes to be created first.
            for n = 1: length(jnodes)
                if isempty(jnodes{n})
                    continue;
                end
                if ~isa(line_nodes{n}, 'Queue')
                    continue;
                end
                % Immediate feedback (self-looping jobs stay in service): either
                % 'all' or a cell array of class indices on the MATLAB side
                if ~isempty(line_nodes{n}.immediateFeedback)
                    if ischar(line_nodes{n}.immediateFeedback)
                        jnodes{n}.setImmediateFeedback(true);
                    elseif iscell(line_nodes{n}.immediateFeedback)
                        for fbc = 1:length(line_nodes{n}.immediateFeedback)
                            jnodes{n}.setImmediateFeedback(jclasses{line_nodes{n}.immediateFeedback{fbc}});
                        end
                    end
                end
                for r = 1:length(line_classes)
                    % Reneging (timer-based patience). Only RENEGING is supported;
                    % BALKING via setPatience is rejected by both MATLAB and the JAR.
                    if ~isempty(line_nodes{n}.patienceDistributions) && ...
                            r <= length(line_nodes{n}.patienceDistributions) && ...
                            ~isempty(line_nodes{n}.patienceDistributions{1, r})
                        patDist = line_nodes{n}.patienceDistributions{1, r};
                        if ~isa(patDist, 'Disabled')
                            impType = ImpatienceType.RENEGING;
                            if ~isempty(line_nodes{n}.impatienceTypes) && ...
                                    r <= length(line_nodes{n}.impatienceTypes) && ...
                                    ~isempty(line_nodes{n}.impatienceTypes{1, r})
                                impType = line_nodes{n}.impatienceTypes{1, r};
                            end
                            jImpType = jline.lang.constant.ImpatienceType.fromID(int32(impType));
                            jnodes{n}.setPatience(jclasses{r}, jImpType, JLINE.from_line_distribution(patDist));
                        end
                    end
                    % Balking (state-based, queue-length/expected-wait thresholds)
                    if ~isempty(line_nodes{n}.balkingStrategies) && ...
                            r <= length(line_nodes{n}.balkingStrategies) && ...
                            ~isempty(line_nodes{n}.balkingStrategies{1, r})
                        balkStrat = line_nodes{n}.balkingStrategies{1, r};
                        switch balkStrat
                            case BalkingStrategy.QUEUE_LENGTH
                                jBalkStrat = jline.lang.constant.BalkingStrategy.QUEUE_LENGTH;
                            case BalkingStrategy.EXPECTED_WAIT
                                jBalkStrat = jline.lang.constant.BalkingStrategy.EXPECTED_WAIT;
                            case BalkingStrategy.COMBINED
                                jBalkStrat = jline.lang.constant.BalkingStrategy.COMBINED;
                            otherwise
                                jBalkStrat = jline.lang.constant.BalkingStrategy.QUEUE_LENGTH;
                        end
                        thresholds = line_nodes{n}.balkingThresholds{1, r};
                        jThresholds = java.util.ArrayList();
                        for ti = 1:length(thresholds)
                            th = thresholds{ti};
                            minJobs = th{1};
                            maxJobs = th{2};
                            if isinf(maxJobs)
                                maxJobs = java.lang.Integer.MAX_VALUE;
                            end
                            jThresholds.add(javaObject('jline.lang.constant.BalkingThreshold', int32(minJobs), int32(maxJobs), th{3}));
                        end
                        jnodes{n}.setBalking(jclasses{r}, jBalkStrat, jThresholds);
                    end
                    % Retrial (orbit + retrial delay distribution)
                    if ~isempty(line_nodes{n}.retrialDelays) && ...
                            r <= length(line_nodes{n}.retrialDelays) && ...
                            ~isempty(line_nodes{n}.retrialDelays{1, r})
                        retDist = line_nodes{n}.retrialDelays{1, r};
                        if ~isa(retDist, 'Disabled')
                            maxAttempts = -1;
                            if ~isempty(line_nodes{n}.retrialMaxAttempts) && r <= length(line_nodes{n}.retrialMaxAttempts)
                                maxAttempts = line_nodes{n}.retrialMaxAttempts(r);
                            end
                            jnodes{n}.setRetrial(jclasses{r}, JLINE.from_line_distribution(retDist), int32(maxAttempts));
                        end
                    end
                    % Orbit impatience (abandonment from the retrial orbit)
                    if ~isempty(line_nodes{n}.orbitImpatienceDistributions) && ...
                            r <= length(line_nodes{n}.orbitImpatienceDistributions) && ...
                            ~isempty(line_nodes{n}.orbitImpatienceDistributions{1, r})
                        orbDist = line_nodes{n}.orbitImpatienceDistributions{1, r};
                        if ~isa(orbDist, 'Disabled')
                            jnodes{n}.setOrbitImpatience(jclasses{r}, JLINE.from_line_distribution(orbDist));
                        end
                    end
                    % Batch rejection probability (retrial queues)
                    if ~isempty(line_nodes{n}.batchRejectProb) && ...
                            r <= length(line_nodes{n}.batchRejectProb) && ...
                            line_nodes{n}.batchRejectProb(r) > 0
                        jnodes{n}.setBatchRejectProbability(jclasses{r}, line_nodes{n}.batchRejectProb(r));
                    end
                end
            end

            % PAS rate function + swap graph transfer; requires classes created first -- see _kb/12-interfaces-and-docs.md
            for n = 1: length(jnodes)
                if isempty(jnodes{n})
                    continue;
                end
                if ~isa(line_nodes{n}, 'Queue') || ...
                        (line_nodes{n}.schedStrategy ~= SchedStrategy.PAS && line_nodes{n}.schedStrategy ~= SchedStrategy.OI)
                    continue;
                end
                if isempty(line_nodes{n}.svcRateFun)
                    line_error(mfilename, sprintf('PAS queue ''%s'' has no service rate function mu(c); set it via setService(@(c) ...).', line_nodes{n}.getName));
                end
                if isinf(line_nodes{n}.cap)
                    line_error(mfilename, sprintf('PAS queue ''%s'' requires a finite capacity for JLINE conversion; set it via setCap(...).', line_nodes{n}.getName));
                end
                jSerFun = JLINE.pas_handle_to_serializablefun(line_nodes{n}.svcRateFun, sn.nclasses, line_nodes{n}.cap);
                jnodes{n}.setServiceRateFunction(jSerFun);
                % Transfer the swap graph if explicitly set (otherwise the JAR
                % defaults to a complete graph at struct refresh, as MATLAB does)
                if ~isempty(line_nodes{n}.swapGraph)
                    jnodes{n}.setSwapGraph(JLINE.from_line_matrix(line_nodes{n}.swapGraph));
                end
            end

            % Transfer heterogeneous server types and their per-(serverType,class)
            % service distributions. Requires classes to be created first.
            for n = 1: length(jnodes)
                if isempty(jnodes{n})
                    continue;
                end
                if ~isa(line_nodes{n}, 'Queue') || isempty(line_nodes{n}.serverTypes)
                    continue;
                end
                % Map MATLAB class names to Java class objects
                % Add each server type with its compatible classes
                jServerTypes = cell(1, length(line_nodes{n}.serverTypes));
                for st = 1:length(line_nodes{n}.serverTypes)
                    serverType = line_nodes{n}.serverTypes{st};
                    jCompat = java.util.ArrayList();
                    compatClasses = serverType.getCompatibleClasses();
                    for cc = 1:length(compatClasses)
                        jCompat.add(jclasses{compatClasses{cc}.index});
                    end
                    jST = javaObject('jline.lang.constant.ServerType', serverType.getName(), int32(serverType.getNumOfServers()), jCompat);
                    jnodes{n}.addServerType(jST);
                    jServerTypes{st} = jST;
                end
                % Set heterogeneous scheduling policy
                switch line_nodes{n}.heteroSchedPolicy
                    case HeteroSchedPolicy.ORDER
                        jHetPol = jline.lang.constant.HeteroSchedPolicy.ORDER;
                    case HeteroSchedPolicy.ALIS
                        jHetPol = jline.lang.constant.HeteroSchedPolicy.ALIS;
                    case HeteroSchedPolicy.ALFS
                        jHetPol = jline.lang.constant.HeteroSchedPolicy.ALFS;
                    case HeteroSchedPolicy.FAIRNESS
                        jHetPol = jline.lang.constant.HeteroSchedPolicy.FAIRNESS;
                    case HeteroSchedPolicy.FSF
                        jHetPol = jline.lang.constant.HeteroSchedPolicy.FSF;
                    case HeteroSchedPolicy.RAIS
                        jHetPol = jline.lang.constant.HeteroSchedPolicy.RAIS;
                    otherwise
                        jHetPol = jline.lang.constant.HeteroSchedPolicy.ORDER;
                end
                jnodes{n}.setHeteroSchedPolicy(jHetPol);
                % Set per-(serverType, class) service distributions
                for st = 1:length(line_nodes{n}.serverTypes)
                    serverType = line_nodes{n}.serverTypes{st};
                    for r = 1:length(line_classes)
                        hetDist = line_nodes{n}.getHeteroService(line_classes{r}, serverType);
                        if ~isempty(hetDist)
                            jnodes{n}.setService(jclasses{r}, jServerTypes{st}, JLINE.from_line_distribution(hetDist));
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
                    % Restore RAND routing for all classes, matching ClosedClass/OpenClass
                    % constructor behavior (initJoinJobClasses sets DISABLED by default)
                    for r = 1 : sn.nclasses
                        jnodes{n}.setRouting(jclasses{r}, jline.lang.constant.RoutingStrategy.RAND);
                    end
                    % Join quorum strategy/count transferred after initJoinJobClasses resets defaults
                    for r = 1 : sn.nclasses
                        if length(line_nodes{n}.input.joinStrategy) >= r && ~isempty(line_nodes{n}.input.joinStrategy{r}) ...
                                && line_nodes{n}.input.joinStrategy{r} == JoinStrategy.PARTIAL
                            jnodes{n}.setStrategy(jclasses{r}, jline.lang.constant.JoinStrategy.PARTIAL);
                        end
                        if length(line_nodes{n}.input.joinRequired) >= r && ~isempty(line_nodes{n}.input.joinRequired{r}) ...
                                && line_nodes{n}.input.joinRequired{r} > 0
                            jnodes{n}.setRequired(jclasses{r}, line_nodes{n}.input.joinRequired{r});
                        end
                    end
                elseif isa(line_nodes{n},"Cache")
                    hitC = line_nodes{n}.server.hitClass;
                    missC = line_nodes{n}.server.missClass;
                    for r = 1 : sn.nclasses
                        % Per-class cache setup gated on what the class actually has, not hitClass alone -- see _kb/09-ldes-and-cache.md
                        hasPop = numel(line_nodes{n}.popularity) >= r ...
                            && ~isempty(line_nodes{n}.popularity{r}) ...
                            && ~isa(line_nodes{n}.popularity{r},'Disabled');
                        hasHit = length(hitC) >= r && full(hitC(r)) > 0;
                        hasMiss = length(missC) >= r && full(missC(r)) > 0;
                        if hasPop
                            jnodes{n}.setRead(jclasses{r}, JLINE.from_line_distribution(line_nodes{n}.popularity{r}));
                        end
                        if hasHit
                            jnodes{n}.setHitClass(jclasses{r}, jclasses{full(hitC(r))});
                        end
                        if hasMiss
                            jnodes{n}.setMissClass(jclasses{r}, jclasses{full(missC(r))});
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
                    % attachRetrievalSystem restores delayed-hit bookkeeping the JAR solvers detect -- see _kb/09-ldes-and-cache.md
                    if ~isempty(line_nodes{n}.retrievalSystemQueueIndices) ...
                            && numEntries(line_nodes{n}.retrievalSystemQueueIndices) > 0
                        rsqi = line_nodes{n}.retrievalSystemQueueIndices;
                        rclasses = line_nodes{n}.server.retrievalClasses; % [nItems x nclasses], 1-based or -1
                        nItemsRS = size(rclasses, 1);
                        keysRS = keys(rsqi);
                        for kk = 1:numel(keysRS)
                            jobinIdx0 = double(keysRS(kk));            % 0-based arrival class index
                            jobinClassObj = jclasses{jobinIdx0 + 1};
                            queueIdxs = rsqi{keysRS(kk)};              % 1-based MATLAB node indices
                            qList = javaObject('java.util.ArrayList');
                            for q = 1:numel(queueIdxs)
                                qList.add(java.lang.Integer(int32(queueIdxs(q) - 1))); % 0-based node index
                            end
                            rcArr = javaArray('jline.lang.JobClass', nItemsRS);
                            for i = 1:nItemsRS
                                rIdx = rclasses(i, jobinIdx0 + 1);
                                if rIdx >= 1
                                    rcArr(i) = jclasses{rIdx};
                                end
                            end
                            jnodes{n}.attachRetrievalSystem(jobinClassObj, qList, rcArr);
                        end
                    end
                elseif isa(line_nodes{n}, "Place") && line_nodes{n}.isQueueing()
                    % QPN embedded queue marshalled after classes exist; setService flips the Place to queueing -- see _kb/12-interfaces-and-docs.md
                    for r = 1:sn.nclasses
                        if numel(line_nodes{n}.serviceProcess) >= r && ~isempty(line_nodes{n}.serviceProcess{r})
                            jdist = JLINE.from_line_distribution(line_nodes{n}.serviceProcess{r});
                            jnodes{n}.setService(jclasses{r}, jdist);
                        end
                    end
                    nsrv = line_nodes{n}.numberOfServers;
                    if ~isinf(nsrv)
                        jnodes{n}.setNumberOfServers(int32(nsrv));
                    end
                    for r = 1:sn.nclasses
                        if numel(line_nodes{n}.departureDiscipline) >= r && line_nodes{n}.departureDiscipline(r) > 0
                            jnodes{n}.setDepartureDiscipline(jclasses{r}, ...
                                jline.lang.constant.DepartureDiscipline.fromID(line_nodes{n}.departureDiscipline(r)));
                        end
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
                        nsrv = line_nodes{n}.numberOfServers(m);
                        if isinf(nsrv)
                            jnodes{n}.setNumberOfServers(jmode, java.lang.Integer(intmax('int32')));
                        else
                            jnodes{n}.setNumberOfServers(jmode, java.lang.Integer(int32(nsrv)));
                        end
                    end
                    % Now set enabling conditions, inhibiting conditions, and firing outcomes
                    jmodes = jnodes{n}.getModes();
                    for m = 1:line_nodes{n}.getNumberOfModes()
                        jmode = jmodes.get(m-1);
                        enabCond = line_nodes{n}.enablingConditions{m};
                        inhibCond = line_nodes{n}.inhibitingConditions{m};
                        firingOut = line_nodes{n}.firingOutcomes{m};
                        % Condition matrices sized at addMode time; bound access by their own row count to avoid over-indexing
                        for r = 1:sn.nclasses
                            for i = 1:length(line_nodes)
                                % Set enabling conditions
                                if i <= size(enabCond, 1) && enabCond(i, r) > 0 && isa(line_nodes{i}, 'Place')
                                    jnodes{n}.setEnablingConditions(jmode, jclasses{r}, jnodes{i}, enabCond(i, r));
                                end
                                % Set inhibiting conditions
                                if i <= size(inhibCond, 1) && inhibCond(i, r) < Inf && isa(line_nodes{i}, 'Place')
                                    jnodes{n}.setInhibitingConditions(jmode, jclasses{r}, jnodes{i}, inhibCond(i, r));
                                end
                                % Set firing outcomes
                                if i <= size(firingOut, 1) && firingOut(i, r) ~= 0
                                    jnodes{n}.setFiringOutcome(jmode, jclasses{r}, jnodes{i}, firingOut(i, r));
                                end
                            end
                        end
                    end
                end
            end

            % Node states transferred before from_line_links to precede the JAR's initDefault validation -- see _kb/12-interfaces-and-docs.md
            for n = 1: length(line_nodes)
                if isempty(jnodes{n})
                    continue;
                end
                if line_nodes{n}.isStateful
                    jnodes{n}.setState(JLINE.from_line_matrix(line_nodes{n}.getState));
                    jnodes{n}.setStateSpace(JLINE.from_line_matrix(line_nodes{n}.getStateSpace));
                    jnodes{n}.setStatePrior(JLINE.from_line_matrix(line_nodes{n}.getStatePrior));
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
                    % Global memory budget passed unrounded (fractional classSize footprints are legitimate)
                    if fcr.globalMaxMemory > 0 && ~isinf(fcr.globalMaxMemory)
                        jfcr.setGlobalMaxMemory(fcr.globalMaxMemory);
                    end
                    % Set per-class max jobs, memory, footprint, weight, drop rules
                    for r = 1:length(line_classes)
                        if length(fcr.classMaxJobs) >= r && fcr.classMaxJobs(r) > 0 && ~isinf(fcr.classMaxJobs(r))
                            jfcr.setClassMaxJobs(jclasses{r}, fcr.classMaxJobs(r));
                        end
                        if length(fcr.classMaxMemory) >= r && fcr.classMaxMemory(r) > 0 && ~isinf(fcr.classMaxMemory(r))
                            jfcr.setClassMaxMemory(jclasses{r}, round(fcr.classMaxMemory(r)));
                        end
                        % classSize passed unrounded, matching the JMT XML writer; default-valued (1) entries skipped as redundant
                        if length(fcr.classSize) >= r && isfinite(fcr.classSize(r)) && fcr.classSize(r) >= 0 && fcr.classSize(r) ~= 1
                            jfcr.setClassSize(jclasses{r}, fcr.classSize(r));
                        end
                        if length(fcr.classWeight) >= r && isfinite(fcr.classWeight(r)) && fcr.classWeight(r) > 0 && fcr.classWeight(r) ~= 1
                            jfcr.setClassWeight(jclasses{r}, fcr.classWeight(r));
                        end
                        if length(fcr.dropRule) >= r
                            % Convert MATLAB DropStrategy numeric to Java DropStrategy enum
                            jDropStrategy = jline.lang.constant.DropStrategy.fromID(fcr.dropRule(r));
                            jfcr.setDropRule(jclasses{r}, jDropStrategy);
                        end
                    end
                    % Transfer linear constraints if set
                    if fcr.hasLinearConstraints()
                        jA = JLINE.from_line_matrix(fcr.constraintA);
                        jb = JLINE.from_line_matrix(fcr.constraintB);
                        jfcr.setLinearConstraints(jA, jb);
                    end
                end
            end

            % CTMC reward handles marshalled as a tabulated TabulatedRewardFunction -- see _kb/12-interfaces-and-docs.md
            if (isstruct(sn) && isfield(sn, 'reward') || isprop(sn, 'reward')) && ~isempty(sn.reward)
                for ri = 1:length(sn.reward)
                    jnetwork.setReward(sn.reward{ri}.name, JLINE.reward_handle_to_tabulatedfun(sn.reward{ri}.fn, sn));
                end
            end

            % Node states transferred before initDefault so it only fills nodes still lacking one
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
            jnetwork.initDefault;
            % Force struct refresh so sn.state reflects updated node states
            jnetwork.setHasStruct(false);

        end

        function jnetwork = line_to_jline(model)
            jnetwork = LINE2JLINE(model);
        end

        function tf = jline_uses_links(jnetwork)
            % True when jline_to_line restores the routing with link(rtorig).
            % State-dependent routing (RROBIN, WRROBIN, JSQ, SQ) cannot go
            % through link(P), which overwrites every strategy with PROB, and a
            % model with no rtorig was never linked in the first place.
            network_nodes = jnetwork.getNodes;
            for n = 1 : network_nodes.size
                output_strategies = network_nodes.get(n-1).getOutputStrategies();
                for m = 1 : output_strategies.size()
                    rs = char(output_strategies.get(m-1).getRoutingStrategy);
                    if any(strcmp(rs, {'RROBIN','WRROBIN','JSQ','SQ'}))
                        tf = false;
                        return;
                    end
                end
            end
            tf = ~isempty(jnetwork.getStruct.rtorig);
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

            % An auto-added ClassSwitch is a PRODUCT of link(P), which re-creates
            % it, so creating it here too collides on the name -- see _kb/12.
            useLinks = JLINE.jline_uses_links(jnetwork);
            autoCS = false(1, network_nodes.size);
            if useLinks
                for n = 1 : network_nodes.size
                    jn = network_nodes.get(n-1);
                    autoCS(n) = isa(jn, 'jline.lang.nodes.ClassSwitch') && jn.autoAdded;
                end
            end

            for n = 1 : network_nodes.size
                if ~isa(network_nodes.get(n-1), 'jline.lang.nodes.ClassSwitch')
                    line_nodes{n} = JLINE.from_jline_node(network_nodes.get(n-1), model, job_classes);
                end
            end

            for n = 1 : job_classes.size
                line_classes{n} = JLINE.from_jline_class(job_classes.get(n-1), model);
            end

            % Deferred signal target association (forJobClass), once all
            % classes exist
            for n = 1 : job_classes.size
                jc = job_classes.get(n-1);
                if isa(jc, 'jline.lang.ClosedSignal') || isa(jc, 'jline.lang.OpenSignal') || isa(jc, 'jline.lang.Signal')
                    jtarget = jc.getTargetJobClass();
                    if ~isempty(jtarget)
                        tname = char(jtarget.getName);
                        for m = 1 : job_classes.size
                            if strcmp(line_classes{m}.name, tname)
                                line_classes{n}.forJobClass(line_classes{m});
                                break;
                            end
                        end
                    end
                end
            end

            for n = 1 : network_nodes.size
                if isa(network_nodes.get(n-1), 'jline.lang.nodes.ClassSwitch') && ~autoCS(n)
                    line_nodes{n} = JLINE.from_jline_node(network_nodes.get(n-1), model, job_classes);
                end
            end

            % Deferred Fork/Join linking: set joinOf on Join nodes
            for n = 1 : network_nodes.size
                jnode = network_nodes.get(n-1);
                if isa(jnode, 'jline.lang.nodes.Join') && ~isempty(jnode.joinOf)
                    forkName = char(jnode.joinOf.getName);
                    for m = 1 : network_nodes.size
                        if ~isempty(line_nodes{m}) && isa(line_nodes{m}, 'Fork') && strcmp(line_nodes{m}.name, forkName)
                            line_nodes{n}.joinOf = line_nodes{m};
                            break;
                        end
                    end
                end
            end

            % Deferred Cache setup: set hitClass, missClass, popularity
            for n = 1 : network_nodes.size
                jnode = network_nodes.get(n-1);
                if isa(jnode, 'jline.lang.nodes.Cache') && isa(line_nodes{n}, 'Cache')
                    cacheNode = line_nodes{n};
                    % hitClass and missClass are stored as index vectors
                    hitClassVec = JLINE.from_jline_matrix(jnode.getHitClass());
                    missClassVec = JLINE.from_jline_matrix(jnode.getMissClass());
                    for r = 1:job_classes.size
                        hitIdx = hitClassVec(r);
                        if hitIdx >= 0 && (hitIdx + 1) <= job_classes.size
                            cacheNode.setHitClass(line_classes{r}, line_classes{hitIdx + 1});
                        end
                        missIdx = missClassVec(r);
                        if missIdx >= 0 && (missIdx + 1) <= job_classes.size
                            cacheNode.setMissClass(line_classes{r}, line_classes{missIdx + 1});
                        end
                    end
                    % Popularity distributions
                    for r = 1:job_classes.size
                        try
                            popDist = jnode.popularityGet(0, r-1);
                            if ~isempty(popDist) && popDist.isDiscrete()
                                matlabDist = JLINE.from_jline_distribution(popDist);
                                if ~isempty(matlabDist)
                                    cacheNode.setRead(line_classes{r}, matlabDist);
                                end
                            end
                        catch
                            % No popularity for this class
                        end
                    end
                end
            end

            for n = 1 : network_nodes.size
                if isempty(line_nodes{n}); continue; end
                JLINE.set_line_service(network_nodes.get(n-1), line_nodes{n}, job_classes, line_classes);
            end

            % Configure Transition modes (distributions, enabling/inhibiting/firing, etc.)
            for n = 1 : network_nodes.size
                jnode = network_nodes.get(n-1);
                if isa(jnode, 'jline.lang.nodes.Transition') && isa(line_nodes{n}, 'Transition')
                    tnode = line_nodes{n};
                    jmodes = jnode.getModes();
                    nmodes = jmodes.size();
                    for m = 1 : nmodes
                        jmode = jmodes.get(m-1);
                        modeName = char(jmode.getName());
                        % addMode if not yet present (Transition starts with one default mode)
                        if m > tnode.getNumberOfModes()
                            tnode.addMode(modeName);
                        else
                            tnode.setModeNames(m, modeName);
                        end
                        % Timing strategy
                        ts = jnode.timingStrategies.get(jmode);
                        if ~isempty(ts)
                            tsName = char(ts.name());
                            if strcmp(tsName, 'IMMEDIATE')
                                tnode.setTimingStrategy(m, TimingStrategy.IMMEDIATE);
                            else
                                tnode.setTimingStrategy(m, TimingStrategy.TIMED);
                            end
                        end
                        % Distribution
                        jdist = jnode.getFiringDistribution(jmode);
                        if ~isempty(jdist)
                            matlabDist = JLINE.from_jline_distribution(jdist);
                            if ~isempty(matlabDist)
                                tnode.setDistribution(m, matlabDist);
                            end
                        end
                        % Number of servers
                        numSrv = jnode.getNumberOfModeServers(jmode);
                        if numSrv == intmax('int32') || numSrv == intmax('int64')
                            tnode.setNumberOfServers(m, Inf);
                        else
                            tnode.setNumberOfServers(m, double(numSrv));
                        end
                        % Firing priority and weight (read from matrices by mode index)
                        if jnode.firingPriorities.getNumElements() > (m-1)
                            tnode.setFiringPriorities(m, jnode.firingPriorities.get(m-1));
                        end
                        if jnode.firingWeights.getNumElements() > (m-1)
                            tnode.setFiringWeights(m, jnode.firingWeights.get(m-1));
                        end
                        % Enabling conditions
                        ecMat = jnode.enablingConditions.get(jmode);
                        if ~isempty(ecMat)
                            nrows = ecMat.getNumRows();
                            ncols = ecMat.getNumCols();
                            for ni = 1 : nrows
                                for ci = 1 : ncols
                                    val = ecMat.get(ni-1, ci-1);
                                    if val > 0 && ni <= length(line_nodes) && isa(line_nodes{ni}, 'Place')
                                        tnode.setEnablingConditions(m, ci, line_nodes{ni}, val);
                                    end
                                end
                            end
                        end
                        % Inhibiting conditions
                        icMat = jnode.inhibitingConditions.get(jmode);
                        if ~isempty(icMat)
                            nrows = icMat.getNumRows();
                            ncols = icMat.getNumCols();
                            for ni = 1 : nrows
                                for ci = 1 : ncols
                                    val = icMat.get(ni-1, ci-1);
                                    if isfinite(val) && val > 0 && ni <= length(line_nodes) && isa(line_nodes{ni}, 'Place')
                                        tnode.setInhibitingConditions(m, ci, line_nodes{ni}, val);
                                    end
                                end
                            end
                        end
                        % Firing outcomes
                        foMat = jnode.firingOutcomes.get(jmode);
                        if ~isempty(foMat)
                            nrows = foMat.getNumRows();
                            ncols = foMat.getNumCols();
                            for ni = 1 : nrows
                                for ci = 1 : ncols
                                    val = foMat.get(ni-1, ci-1);
                                    if val ~= 0 && ni <= length(line_nodes)
                                        tnode.setFiringOutcome(m, ci, line_nodes{ni}, val);
                                    end
                                end
                            end
                        end
                    end
                end
            end

            if useLinks
                % Use link() method
                model = JLINE.from_jline_links(model, jnetwork);
            else
                % Do not use link() method
                model = JLINE.from_jline_routing(model, jnetwork);
            end

            % Restore initial state on Place nodes after linking
            for n = 1 : network_nodes.size
                jnode = network_nodes.get(n-1);
                if isa(jnode, 'jline.lang.nodes.Place') && ~isempty(line_nodes{n}) && isa(line_nodes{n}, 'Place')
                    jst = jnode.getState();
                    if ~isempty(jst) && ~jst.isEmpty()
                        line_nodes{n}.setState(JLINE.from_jline_matrix(jst));
                    end
                end
            end

            % FCR transfer, reverse of the from_line_network block above
            jregions = jnetwork.getRegions();
            for f = 1 : jregions.size()
                jfcr = jregions.get(f-1);
                jrnodes = jfcr.getNodes();
                regionNodes = cell(1, jrnodes.size());
                for i = 1 : jrnodes.size()
                    rname = char(jrnodes.get(i-1).getName());
                    for n = 1 : length(line_nodes)
                        if ~isempty(line_nodes{n}) && strcmp(line_nodes{n}.getName(), rname)
                            regionNodes{i} = line_nodes{n};
                            break;
                        end
                    end
                end
                fcr = model.addRegion(regionNodes);
                gmj = jfcr.getGlobalMaxJobs();
                if gmj > 0
                    fcr.setGlobalMaxJobs(gmj);
                end
                gmm = jfcr.getGlobalMaxMemory();
                if gmm > 0
                    fcr.setGlobalMaxMemory(gmm);
                end
                for r = 1 : length(line_classes)
                    jclass = job_classes.get(r-1);
                    cmj = jfcr.getClassMaxJobs(jclass);
                    if cmj > 0
                        fcr.setClassMaxJobs(line_classes{r}, cmj);
                    end
                    cmm = jfcr.getClassMaxMemory(jclass);
                    if cmm > 0
                        fcr.setClassMaxMemory(line_classes{r}, cmm);
                    end
                    cw = jfcr.getClassWeight(jclass);
                    if isfinite(cw) && cw > 0 && cw ~= 1
                        fcr.setClassWeight(line_classes{r}, cw);
                    end
                    cs = jfcr.getClassSize(jclass);
                    if isfinite(cs) && cs >= 0 && cs ~= 1
                        fcr.setClassSize(line_classes{r}, cs);
                    end
                    % Java getID() uses the MATLAB numeric ids (fromID accepts
                    % them on the forward path), so the id round-trips directly.
                    fcr.setDropRule(line_classes{r}, jfcr.getDropStrategy(jclass).getID());
                end
                jlin = jfcr.getLinearConstraints();
                if ~isempty(jlin)
                    fcr.setConstraint(JLINE.from_jline_matrix(jlin(1)), JLINE.from_jline_matrix(jlin(2)));
                end
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

        function out = from_jline_matrix_list(jlist)
            % java.util.List<Matrix> -> 1-by-n cell of double matrices
            out = {};
            if isempty(jlist), return; end
            n = double(jlist.size());
            out = cell(1, n);
            for k = 1:n
                out{k} = JLINE.from_jline_matrix(jlist.get(k-1));
            end
        end

        function out = from_jline_matrixcell(jcell)
            % jline.util.matrix.MatrixCell -> 1-by-n cell of double matrices
            out = {};
            if isempty(jcell), return; end
            n = double(jcell.size());
            out = cell(1, n);
            for k = 1:n
                out{k} = JLINE.from_jline_matrix(jcell.get(k-1));
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
                        % No band around 1e8/3 here: the JAR carries the same
                        % Immediate constant (1/FineTol) as MATLAB, and a rate
                        % of 33333333.33 is a legitimate Exp(3e-8), e.g. three
                        % chained immediate demands in an LN layer.
                        if (val >= 2147483647 - 1) % Integer.MAX_VALUE with -1 tolerance
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

        function jdist = from_line_lqn_dist(ptype, dmean, dscv, dparams, dproc)
            % JDIST = FROM_LINE_LQN_DIST(PTYPE, DMEAN, DSCV, DPROC)
            % Rebuild a JAR distribution from the LayeredNetworkStruct fields of
            % an LQN distribution, keeping its variability. Passing the mean
            % alone to a setSomething(double) overload rebuilds it as Exp(1/mean)
            % with SCV 1, which silently discards everything above the first
            % moment: an Erlang setup and an exponential one of the same mean
            % would reach the JAR as the same process.
            switch ptype
                case ProcessType.IMMEDIATE
                    jdist = jline.lang.processes.Immediate;
                case ProcessType.EXP
                    jdist = jline.lang.processes.Exp(1/dmean);
                case ProcessType.ERLANG
                    jdist = jline.lang.processes.Erlang.fitMeanAndSCV(dmean, dscv);
                case ProcessType.HYPEREXP
                    % Rebuilt from its own (p,lambda1,lambda2), not a two-moment refit -- see _kb/12-interfaces-and-docs.md
                    if ~isempty(dparams) && length(dparams) >= 3
                        jdist = jline.lang.processes.HyperExp(dparams(1), dparams(2), dparams(3));
                    else
                        jdist = jline.lang.processes.HyperExp.fitMeanAndSCV(dmean, dscv);
                    end
                case ProcessType.COXIAN
                    jdist = jline.lang.processes.Coxian.fitMeanAndSCV(dmean, dscv);
                case ProcessType.APH
                    jdist = jline.lang.processes.APH.fitMeanAndSCV(dmean, dscv);
                case {ProcessType.PH, ProcessType.MAP}
                    if ~isempty(dproc)
                        if ptype == ProcessType.PH
                            jdist = jline.lang.processes.PH(JLINE.from_line_matrix(dproc{1}), JLINE.from_line_matrix(dproc{2}));
                        else
                            jdist = jline.lang.processes.MAP(JLINE.from_line_matrix(dproc{1}), JLINE.from_line_matrix(dproc{2}));
                        end
                    else
                        jdist = jline.lang.processes.Exp(1/dmean);
                    end
                case ProcessType.DET
                    jdist = jline.lang.processes.Det(dmean);
                otherwise
                    % Any other type is carried by mean and SCV, which is exact
                    % for SCV 1 and a two-moment fit otherwise.
                    if abs(dscv-1) < GlobalConstants.FineTol
                        jdist = jline.lang.processes.Exp(1/dmean);
                    else
                        jdist = jline.lang.processes.APH.fitMeanAndSCV(dmean, dscv);
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
            % JLINE indexes LayeredNetworkStruct elements from 0 and carries no
            % padding row or column; MATLAB indexes them from 1. This is the only
            % seam between the two conventions, so every element index crossing it
            % is shifted by one here: map keys are looked up at key-1, and index
            % VALUES (parent, callpair, tasksof/entriesof/actsof/callsof) come back
            % +1. Java's -1 "unset" parent becomes MATLAB's 0.
            % See _kb/04-networkstruct.md and _kb/07-cross-language-parity.md.
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
                lsn.tasksof{h,1} = JLINE.shift_idx(JLINE.arraylist_to_matrix(jlsn.tasksof.get(uint32(h-1))))';
            end
            for t=1:jlsn.ntasks
                lsn.entriesof{lsn.tshift+t,1} = JLINE.shift_idx(JLINE.arraylist_to_matrix(jlsn.entriesof.get(uint32(jlsn.tshift+t-1))))';
            end
            for t=1:(jlsn.ntasks+jlsn.nentries)
                lsn.actsof{lsn.tshift+t,1} = JLINE.shift_idx(JLINE.arraylist_to_matrix(jlsn.actsof.get(uint32(jlsn.tshift+t-1))))';
            end
            for a=1:jlsn.nacts
                lsn.callsof{lsn.ashift+a,1} = JLINE.shift_idx(JLINE.arraylist_to_matrix(jlsn.callsof.get(uint32(jlsn.ashift+a-1))))';
            end
            for i = 1:jlsn.sched.size
                % A class's Constant properties do NOT accept dynamic field
                % access: SchedStrategy.(name) resolves the identifier
                % 'SchedStrategy.' as a class and errors whatever name holds.
                % Go through the same .name.toCharArray' + fromText pair the
                % flat sn.sched conversion below uses.
                jsched = jlsn.sched.get(uint32(i-1));
                lsn.sched(i,1) = SchedStrategy.fromText(jsched.name.toCharArray');
            end
            for i = 1:jlsn.names.size
                lsn.names{i,1} = jlsn.names.get(uint32(i-1));
                lsn.hashnames{i,1} = jlsn.hashnames.get(uint32(i-1));
            end
            lsn.mult = JLINE.from_jline_matrix(jlsn.mult);
            lsn.mult = lsn.mult(1:lsn.eshift)';
            lsn.maxmult = JLINE.from_jline_matrix(jlsn.maxmult);
            lsn.maxmult = lsn.maxmult(1:lsn.eshift)';

            lsn.repl = JLINE.from_jline_matrix(jlsn.repl)';
            lsn.type = JLINE.from_jline_matrix(jlsn.type)';
            % parent holds element indices: shift them, and map Java's -1 to 0
            lsn.parent = JLINE.shift_idx(JLINE.from_jline_matrix(jlsn.parent));
            lsn.nitems = JLINE.from_jline_matrix(jlsn.nitems);
            % Ensure proper column vector format matching MATLAB's (nhosts+ntasks+nentries) x 1
            if isrow(lsn.nitems)
                lsn.nitems = lsn.nitems';
            end
            % Ensure correct size
            expectedSize = lsn.nhosts + lsn.ntasks + lsn.nentries;
            if length(lsn.nitems) < expectedSize
                lsn.nitems(expectedSize,1) = 0;
            elseif length(lsn.nitems) > expectedSize
                lsn.nitems = lsn.nitems(1:expectedSize);
            end
            lsn.replacestrat = JLINE.from_jline_matrix(jlsn.replacestrat)';
            for i = 1:jlsn.callnames.size
                lsn.callnames{i,1} = jlsn.callnames.get(uint32(i-1));
                lsn.callhashnames{i,1} = jlsn.callhashnames.get(uint32(i-1));
            end
            for i = 1:jlsn.calltype.size % calltype may be made into a matrix in Java
                % CallType.(ct) is the same unsupported dynamic access as
                % SchedStrategy.(...) above; switch on the enum's own name
                ct = jlsn.calltype.get(uint32(i-1)).name.toCharArray';
                switch ct
                    case 'SYNC'
                        lsn.calltype(i) = CallType.SYNC;
                    case 'ASYNC'
                        lsn.calltype(i) = CallType.ASYNC;
                    case 'FWD'
                        lsn.calltype(i) = CallType.FWD;
                    otherwise
                        line_error(mfilename, sprintf('Unknown call type %s in the Java LayeredNetworkStruct.', ct));
                end
            end
            lsn.calltype = sparse(lsn.calltype');
            % callpair rows are calls and its entries are element indices
            lsn.callpair = JLINE.shift_idx(JLINE.from_jline_matrix(jlsn.callpair));
            if isempty(lsn.callpair)
                lsn.callpair=[];
            end
            lsn.actpretype = sparse(JLINE.from_jline_matrix(jlsn.actpretype)');
            lsn.actposttype = sparse(JLINE.from_jline_matrix(jlsn.actposttype)');
            lsn.graph = JLINE.from_jline_matrix(jlsn.graph);
            lsn.dag = JLINE.from_jline_matrix(jlsn.dag);
            lsn.taskgraph = sparse(JLINE.from_jline_matrix(jlsn.taskgraph));
            lsn.replygraph = logical(JLINE.from_jline_matrix(jlsn.replygraph));
            lsn.iscache = JLINE.from_jline_matrix(jlsn.iscache);
            % Ensure proper column vector format matching MATLAB's (nhosts+ntasks) x 1
            % The JAR row is 1 x nidx indexed by the GLOBAL element index, and
            % since 778978b66 that index is 0-BASED: element 1 of the row is
            % host 1, not a pad. Dropping it, as the 1-based era required,
            % shifted every host/task left by one and reported the cache task
            % one slot early (lsnDebug: LINE [0 0 0 1] vs JLINE [0 0 1 0] on
            % lcq_singlehost). Hosts and tasks occupy the FIRST nhosts+ntasks
            % columns, so those are the ones to take.
            expectedCacheSize = lsn.nhosts + lsn.ntasks;
            if isrow(lsn.iscache)
                if length(lsn.iscache) > expectedCacheSize
                    lsn.iscache = lsn.iscache(1:expectedCacheSize)'; % drop the trailing element columns
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
            lsn.iscaller = full(JLINE.from_jline_matrix(jlsn.iscaller));
            lsn.issynccaller = full(JLINE.from_jline_matrix(jlsn.issynccaller));
            lsn.isasynccaller = full(JLINE.from_jline_matrix(jlsn.isasynccaller));
            lsn.isref = JLINE.from_jline_matrix(jlsn.isref)';
        end

        function out = shift_idx(idxs)
            % Convert JLINE 0-based element (or call) indices to MATLAB 1-based
            % ones. Java marks "unset" with -1, which MATLAB spells as 0, so the
            % same +1 carries both. An empty input stays empty.
            out = idxs;
            if isempty(out)
                return;
            end
            out = out + 1;
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
                stationFunMap = configureDictionary('string','cell');
                while entryIter.hasNext()
                    entry = entryIter.next();
                    stationName = char(entry.getKey().getName());
                    try
                        jfun = entry.getValue();
                        if ~isempty(jfun)
                            stationFunMap{stationName} = jfun;
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
                        jfun = stationFunMap{stationName};
                        % Create a MATLAB function handle that calls the Java apply() method
                        sn.cdscaling{i} = @(ni) JLINE.call_java_cdscaling(jfun, ni);
                    else
                        sn.cdscaling{i} = @(ni) 1;
                    end
                end
            else
                sn.cdscaling = cell(sn.nstations, 0);
            end

            % joint-dependence handles eta_i(n) (non-product-form), twin of the
            % cdscaling readback above.
            if isprop(jsn, 'jdscaling') && ~isempty(jsn.jdscaling) && jsn.jdscaling.size() > 0
                sn.jdscaling = cell(sn.nstations, 1);
                entrySet = jsn.jdscaling.entrySet();
                entryIter = entrySet.iterator();
                stationFunMap = configureDictionary('string','cell');
                while entryIter.hasNext()
                    entry = entryIter.next();
                    stationName = char(entry.getKey().getName());
                    try
                        jfun = entry.getValue();
                        if ~isempty(jfun)
                            stationFunMap{stationName} = jfun;
                        end
                    catch
                    end
                end
                for i = 1:sn.nstations
                    jstation = jstations.get(i-1);
                    stationName = char(jstation.getName());
                    if isKey(stationFunMap, stationName)
                        jfun = stationFunMap{stationName};
                        sn.jdscaling{i} = @(ni) JLINE.call_java_cdscaling(jfun, ni);
                    else
                        sn.jdscaling{i} = @(ni) 1;
                    end
                end
            else
                sn.jdscaling = cell(sn.nstations, 0);
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
                            case 'SQ'
                                sn.routing(i,j) = RoutingStrategy.SQ;
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
                            % mode-indexed numeric ProcessType ids, as refreshPetriNetNodes
                            % builds it; Mode has no toString, so it cannot key a map
                            modeIdx = [];
                            procIds = [];
                            keys = jparam.firingprocid.keySet.iterator;
                            while keys.hasNext
                                key = keys.next;
                                proc = jparam.firingprocid.get(key);
                                modeIdx(end+1) = double(key.getIndex()); %#ok<AGROW>
                                procIds(end+1) = ProcessType.fromText(char( ...
                                    jline.lang.constant.ProcessType.toText(proc))); %#ok<AGROW>
                            end
                            fpid = -ones(1, max([modeIdx, 0]));
                            fpid(modeIdx) = procIds;
                            sn.nodeparam{i}.firingprocid = fpid;
                        end
                        if ~isempty(jparam.firingphases)
                            sn.nodeparam{i}.firingphases = JLINE.from_jline_matrix(jparam.firingphases);
                        end
                        if ~isempty(jparam.fireweight)
                            sn.nodeparam{i}.fireweight = JLINE.from_jline_matrix(jparam.fireweight);
                        end
                        nmodes = double(jparam.nmodes);
                        if nmodes == 0 && ~isempty(jparam.modenames)
                            nmodes = double(jparam.modenames.size());
                        end
                        sn.nodeparam{i}.nmodes = nmodes;
                        if ~isempty(jparam.modenames)
                            mn = cell(1, double(jparam.modenames.size()));
                            for m = 1:numel(mn)
                                mn{m} = char(jparam.modenames.get(m-1));
                            end
                            sn.nodeparam{i}.modenames = mn;
                        end
                        sn.nodeparam{i}.enabling = JLINE.from_jline_matrix_list(jparam.enabling);
                        sn.nodeparam{i}.inhibiting = JLINE.from_jline_matrix_list(jparam.inhibiting);
                        sn.nodeparam{i}.firing = JLINE.from_jline_matrix_list(jparam.firing);
                        if ~isempty(jparam.nmodeservers)
                            sn.nodeparam{i}.nmodeservers = JLINE.from_jline_matrix(jparam.nmodeservers);
                        end
                        if ~isempty(jparam.firingprio)
                            sn.nodeparam{i}.firingprio = JLINE.from_jline_matrix(jparam.firingprio);
                        end
                        if ~isempty(jparam.timing)
                            tm = zeros(1, double(jparam.timing.size()));
                            for m = 1:numel(tm)
                                if strcmp(char(jparam.timing.get(m-1).toString), 'IMMEDIATE')
                                    tm(m) = TimingStrategy.IMMEDIATE;
                                else
                                    tm(m) = TimingStrategy.TIMED;
                                end
                            end
                            sn.nodeparam{i}.timing = tm;
                        end
                        % mode-keyed maps: Mode.getIndex is the only stable key
                        if ~isempty(jparam.firingproc)
                            fproc = cell(1, nmodes);
                            it = jparam.firingproc.keySet.iterator;
                            while it.hasNext
                                key = it.next;
                                fproc{double(key.getIndex())} = ...
                                    JLINE.from_jline_matrixcell(jparam.firingproc.get(key));
                            end
                            sn.nodeparam{i}.firingproc = fproc;
                        end
                        if ~isempty(jparam.firingpie)
                            fpie = cell(1, nmodes);
                            it = jparam.firingpie.keySet.iterator;
                            while it.hasNext
                                key = it.next;
                                fpie{double(key.getIndex())} = ...
                                    JLINE.from_jline_matrix(jparam.firingpie.get(key));
                            end
                            sn.nodeparam{i}.firingpie = fpie;
                        end
                        % g_m(marking): the Java lambda cannot become a MATLAB
                        % handle, so wrap it and marshal the marking per call
                        fdep = cell(1, nmodes);
                        if ~isempty(jparam.firingdep)
                            for m = 1:min(nmodes, double(jparam.firingdep.size()))
                                jf = jparam.firingdep.get(m-1);
                                if ~isempty(jf)
                                    fdep{m} = @(M) double(jf.apply(JLINE.from_line_matrix(M)));
                                end
                            end
                        end
                        sn.nodeparam{i}.firingdep = fdep;
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

            if ~isempty(jsn.sync)
                jsync = jsn.sync;
                sn.sync = cell(jsync.size, 1);
                for i = 1:jsync.size
                    jsync_i = jsync.get(uint32(i-1));
                    sn.sync{i,1} = struct('active',cell(1),'passive',cell(1));

                    jactive = jsync_i.active.get(uint32(0));
                    jpassive = jsync_i.passive.get(uint32(0));

                    % Assumes prob is a value, not a Java lambda function.
                    % The mapping is by NAME, once, for every member of the
                    % Java enum: the two hand-written switches this replaces
                    % listed only INIT/LOCAL/ARV/DEP/PHASE/READ/STAGE and left
                    % sn.sync{i}.active{1} UNASSIGNED for anything else, so a
                    % reneging, retrial, polling, breakdown or (now) tagged
                    % model crossed the bridge with a hole in its sync list.
                    sn.sync{i,1}.active{1} = JLINE.eventFromJava(jactive);
                    sn.sync{i,1}.passive{1} = JLINE.eventFromJava(jpassive);
                end
            else
                sn.sync = {};
            end
        end

        function ev = eventFromJava(jev)
            % EV = EVENTFROMJAVA(JEV) build the MATLAB Event of a jline.lang.Event.
            %
            % The Java and MATLAB EventType numberings disagree (Java ordinals
            % start at INIT = 0, MATLAB at INIT = -1), so the two are matched by
            % NAME, as the rest of the codebase does. Every member of the Java
            % enum is listed: an unmapped one must raise rather than leave the
            % synchronization silently unassigned.
            nm = jev.getEvent.name.toCharArray';
            switch nm
                case 'INIT',    et = EventType.INIT;
                case 'LOCAL',   et = EventType.LOCAL;
                case 'ARV',     et = EventType.ARV;
                case 'DEP',     et = EventType.DEP;
                case 'PHASE',   et = EventType.PHASE;
                case 'READ',    et = EventType.READ;
                case 'STAGE',   et = EventType.STAGE;
                case 'ENABLE',  et = EventType.ENABLE;
                case 'FIRE',    et = EventType.FIRE;
                case 'PRE',     et = EventType.PRE;
                case 'POST',    et = EventType.POST;
                case 'RENEGE',  et = EventType.RENEGE;
                case 'RETRY',   et = EventType.RETRY;
                case 'SWITCH',  et = EventType.SWITCH;
                case 'FAILURE', et = EventType.FAILURE;
                case 'REPAIR',  et = EventType.REPAIR;
                case 'START',   et = EventType.START;
                case 'PREEMPT', et = EventType.PREEMPT;
                otherwise
                    line_error(mfilename, sprintf(['The JAR declares event type ''%s'', which this bridge cannot map ' ...
                        'to a MATLAB EventType. Add it to JLINE.eventFromJava rather than letting the ' ...
                        'synchronization cross unassigned.'], nm));
            end
            ev = Event(et, jev.getNode+1, jev.getJobClass+1, ...
                jev.getProb, JLINE.from_jline_matrix(jev.getState), ...
                jev.getT, jev.getJob);
        end

        function [QN,UN,RN,WN,AN,TN] = arrayListToResults(alist)
            switch class(alist)
                case 'jline.solvers.LayeredNetworkAvgTable'
                    QN = JLINE.arraylist_to_matrix(alist.getQLen());
                    UN = JLINE.arraylist_to_matrix(alist.getUtil());
                    RN = JLINE.arraylist_to_matrix(alist.getRespT());
                    WN = JLINE.arraylist_to_matrix(alist.getResidT());
                    AN = JLINE.arraylist_to_matrix(alist.getArvR());
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
                            case 'config'
                                if isfield(options.config,'eventcache')
                                    solverOptions.config.eventcache = options.config.eventcache;
                                end
                                if isfield(options.config,'fork_join')
                                    solverOptions.config.fork_join = options.config.fork_join;
                                end
                                if isfield(options.config,'highvar')
                                    solverOptions.config.highvar = options.config.highvar;
                                end
                                if isfield(options.config,'multiserver')
                                    solverOptions.config.multiserver = options.config.multiserver;
                                end
                                if isfield(options.config,'np_priority')
                                    solverOptions.config.np_priority = options.config.np_priority;
                                end
                                if isfield(options.config,'warmupfrac') && ~isempty(options.config.warmupfrac) ...
                                        && options.config.warmupfrac > 0
                                    % SSA warmup discard (mean estimates + CI batch means)
                                    solverOptions.config.warmupfrac = java.lang.Double(options.config.warmupfrac);
                                end
                                % SolverMAM 'bgchain'. This whitelist is the whole
                                % bridge: a config field not named here is dropped
                                % SILENTLY, so lang='java' would answer the default
                                % while reporting the requested method, which is
                                % indistinguishable from the option having no effect.
                                % SolverMAM 'bgchain'. These three have no declared
                                % field on SolverOptions.Config: they ride its
                                % additionalParams map, which is what the JAR reads
                                % with config.get(name), so they must be PUT rather
                                % than assigned. The map handle needs a TEMPORARY --
                                % MATLAB parses solverOptions.config.put(...) as
                                % nested field indexing on a Java object and fails
                                % with "Dot indexing is not supported", so binding
                                % the Config to a variable first is what makes the
                                % call a method call.
                                % SolverCTMC transient path: a declared field
                                % each, so they are assigned rather than put.
                                if isfield(options.config,'transient_method') && ~isempty(options.config.transient_method)
                                    solverOptions.config.transient_method = options.config.transient_method;
                                end
                                if isfield(options.config,'fau_epsilon') && ~isempty(options.config.fau_epsilon)
                                    solverOptions.config.fau_epsilon = options.config.fau_epsilon;
                                end
                                if isfield(options.config,'fau_delta') && ~isempty(options.config.fau_delta)
                                    solverOptions.config.fau_delta = options.config.fau_delta;
                                end
                                if isfield(options.config,'fau_ngrid') && ~isempty(options.config.fau_ngrid)
                                    solverOptions.config.fau_ngrid = int32(options.config.fau_ngrid);
                                end
                                jconfig = solverOptions.config;
                                if isfield(options.config,'bgaggr') && ~isempty(options.config.bgaggr)
                                    jconfig.put('bgaggr', java.lang.Integer(int32(options.config.bgaggr)));
                                end
                                if isfield(options.config,'bgstates_max') && ~isempty(options.config.bgstates_max)
                                    jconfig.put('bgstates_max', java.lang.Integer(int32(options.config.bgstates_max)));
                                end
                                if isfield(options.config,'qbdphases_max') && ~isempty(options.config.qbdphases_max)
                                    jconfig.put('qbdphases_max', java.lang.Integer(int32(options.config.qbdphases_max)));
                                end
                                % bgenv='full' is a MATLAB-ONLY oracle (see
                                % mam_bgchain_envfull.m); the JAR has no such key, so
                                % bridging it would silently run the lump instead.
                                if isfield(options.config,'bgenv') && ~isempty(options.config.bgenv) ...
                                        && strcmpi(char(options.config.bgenv),'full')
                                    line_error(mfilename,['options.config.bgenv=''full'' is implemented in ' ...
                                        'MATLAB only (mam_bgchain_envfull.m) and has no counterpart in the JAR. ' ...
                                        'Run it with lang=''matlab'', or drop it: it is an oracle for the lumping ' ...
                                        'in mam_bgchain_env and returns the same numbers on every model.']);
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

        function [qns] = SolverQNS(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.QNS);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            qns = jline.solvers.wrappers.qns.SolverQNS(network_object, solverOptions);
        end

        function [mam] = SolverMAM(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.MAM);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            mam = jline.solvers.mam.SolverMAM(network_object, solverOptions);
        end

        function [ag] = SolverAG(network_object, options)
            % The agent-based (RCAT) solver. Its options carry the truncation
            % level of an open agent and the execution backend, neither of which
            % SolverOptions('MAM') has, so it builds AGOptions rather than the
            % generic container.
            solverOptions = jline.solvers.ag.AGOptions();
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            ag = jline.solvers.ag.SolverAG(network_object, solverOptions);
        end

        function [jmt] = SolverJMT(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.JMT);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            jmt = jline.solvers.wrappers.jmt.SolverJMT(network_object, solverOptions);
        end

        function [ctmc] = SolverCTMC(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.CTMC);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            ctmc = jline.solvers.ctmc.SolverCTMC(network_object,solverOptions);
        end

        function [infGen, eventFilt, syncInfo, stateSpace, nodeStateSpace] = getSymbolicGenerator(ctmc, invertSymbol)
            % [INFGEN, EVENTFILT, SYNCINFO, STATESPACE, NODESTATESPACE] = GETSYMBOLICGENERATOR(CTMC, INVERTSYMBOL)
            % Symbolic infinitesimal generator of a JLINE SolverCTMC object, with
            % each event filtration normalized by its minimum positive rate and
            % scaled by a symbolic variable x1..xE, as in the native
            % SolverCTMC.getSymbolicGenerator. Coefficient matrices are computed
            % by the JAR; symbolic objects are rebuilt with the Symbolic Toolbox.
            if nargin<2
                invertSymbol = false;
            end
            if ~exist('sym')
                line_error(mfilename,'This method requires MATLAB''s Symbolic Toolbox.');
            end
            res = ctmc.getSymbolicGenerator(invertSymbol);
            stateSpace = JLINE.from_jline_matrix(res.stateSpace);
            n = size(stateSpace,1);
            nEvents = res.eventFilt.size();
            infGen = sym(zeros(n));
            eventFilt = cell(1, nEvents);
            for e = 1:nEvents
                symName = res.symbols.get(e-1);
                if ~isempty(symName)
                    Fe = JLINE.from_jline_matrix(res.eventFilt.get(e-1));
                    if invertSymbol
                        eventFilt{e} = Fe / sym(char(symName),'real');
                    else
                        eventFilt{e} = Fe * sym(char(symName),'real');
                    end
                    infGen = infGen + eventFilt{e};
                end
            end
            infGen = ctmc_makeinfgen(infGen);
            syncInfo = res.syncInfo;
            nodeStateSpace = cell(1, res.nodeStateSpace.size());
            for i = 1:res.nodeStateSpace.size()
                nodeStateSpace{i} = JLINE.from_jline_matrix(res.nodeStateSpace.get(i-1));
            end
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

        function [ldes] = SolverLDES(network_object, options)
            % Create LDES-specific options object
            ldesOptions = jline.solvers.ldes.LDESOptions();
            if nargin>1
                % Copy standard options
                ldesOptions.samples = options.samples;
                ldesOptions.seed = options.seed;
                % Java backend silenced at SILENT/STD to avoid duplicating the MATLAB summary lines -- see _kb/12-interfaces-and-docs.md
                if isfield(options, 'verbose')
                    switch options.verbose
                        case {VerboseLevel.DEBUG}
                            ldesOptions.verbose = ldesOptions.verbose.DEBUG;
                        otherwise
                            ldesOptions.verbose = ldesOptions.verbose.SILENT;
                    end
                end
                % Parse confint
                [confintEnabled, confintLevel] = Solver.parseConfInt(options.confint);
                if confintEnabled
                    ldesOptions.confint = confintLevel;
                else
                    ldesOptions.confint = 0;
                end
                % Pass timespan for transient analysis
                if isfield(options, 'timespan') && length(options.timespan) >= 2
                    ldesOptions.timespan = options.timespan;
                end
                % Warm-start initial placement (station-major vector, e.g. set
                % by @SolverLDES/initFromSolver from an auxiliary solver)
                if isfield(options, 'init_sol') && ~isempty(options.init_sol)
                    ldesOptions.init_sol = JLINE.from_line_matrix(options.init_sol);
                end
                % Pass LDES-specific options if configured
                if isfield(options, 'config')
                    % Transient detection options
                    if isfield(options.config, 'tranfilter')
                        ldesOptions.tranfilter = options.config.tranfilter;
                    end
                    if isfield(options.config, 'mserbatch')
                        ldesOptions.mserbatch = options.config.mserbatch;
                    end
                    if isfield(options.config, 'warmupfrac')
                        ldesOptions.warmupfrac = options.config.warmupfrac;
                    end
                    % Confidence interval options
                    if isfield(options.config, 'cimethod')
                        ldesOptions.cimethod = options.config.cimethod;
                    end
                    if isfield(options.config, 'obmoverlap')
                        ldesOptions.obmoverlap = options.config.obmoverlap;
                    end
                    if isfield(options.config, 'ciminbatch')
                        ldesOptions.ciminbatch = options.config.ciminbatch;
                    end
                    if isfield(options.config, 'ciminobs')
                        ldesOptions.ciminobs = options.config.ciminobs;
                    end
                    if isfield(options.config, 'spectralLowFreqFrac')
                        ldesOptions.spectralLowFreqFrac = options.config.spectralLowFreqFrac;
                    end
                    % Convergence options
                    if isfield(options.config, 'cnvgon')
                        ldesOptions.cnvgon = options.config.cnvgon;
                    end
                    if isfield(options.config, 'cnvgtol')
                        ldesOptions.cnvgtol = options.config.cnvgtol;
                    end
                    if isfield(options.config, 'cnvgbatch')
                        ldesOptions.cnvgbatch = options.config.cnvgbatch;
                    end
                    if isfield(options.config, 'cnvgchk')
                        ldesOptions.cnvgchk = options.config.cnvgchk;
                    end
                end
            end
            ldes = jline.solvers.ldes.SolverLDES(network_object, ldesOptions);
        end

        function [mva] = SolverMVA(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.MVA);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            mva = jline.solvers.mva.SolverMVA(network_object, solverOptions);
        end

        function [ba] = SolverBA(network_object, options)
            solverOptions = jline.solvers.ba.SolverBA.defaultOptions();
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            ba = jline.solvers.ba.SolverBA(network_object, solverOptions);
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
                % AUTO carries its selection token in options.method, while the
                % JAR keeps it in selectionMethod; without this it stays default.
                if isfield(options,'method') && ~isempty(options.method)
                    solverOptions.selectionMethod = options.method;
                end
            end
            auto = jline.solvers.auto.SolverAUTO(network_object, solverOptions);
        end

        function streamOpts = StreamingOptions(varargin)
            % STREAMINGOPTIONS Create Java StreamingOptions for SSA/LDES stream() method
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

        function [ln] = SolverLN(layered_network_object, options, layerSolverType)
            % LN = SOLVERLN(LAYERED_NETWORK_OBJECT, OPTIONS, LAYERSOLVERTYPE)
            %
            % LAYERSOLVERTYPE is the JAR SolverType of the LAYER solvers the
            % caller asked for. Without it the JAR falls back to its own
            % DefaultSolverFactory, which is SolverMVA at every layer, so
            % LN(model, @(l)NC(l,...)) under lang='java' silently answered with
            % MVA layers -- the same defect the python bridge carries a fix for
            % in PYLINE.SolverLN. The layer solver is what the fixed point is a
            % fixed point OF, so substituting one answers a different question.
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.LN);
            if nargin>1
                solverOptions = JLINE.parseSolverOptions(solverOptions, options);
            end
            if nargin<3 || isempty(layerSolverType)
                ln = jline.solvers.ln.SolverLN(layered_network_object, solverOptions);
            else
                ln = jline.solvers.ln.SolverLN(layered_network_object, layerSolverType, solverOptions);
            end
        end

        function layerSolverType = lnLayerSolverType(solver)
            % LAYERSOLVERTYPE = LNLAYERSOLVERTYPE(SOLVER)
            % Resolve the JAR SolverType of the layer solver a MATLAB SolverLN
            % was built with, refusing anything the JAR has no layer factory
            % for. Empty means "no factory recorded", i.e. the caller kept the
            % default and the JAR may keep its own. The factory is probed on a
            % throwaway network, the same resolution CPPLINE.lnLayerSolver and
            % PYLINE.lnLayerSolverName perform: under lang='java' no MATLAB
            % layer is constructed, so there is no self.solvers to read it off.
            layerSolverType = [];
            if ~isprop(solver, 'solverFactory') || isempty(solver.solverFactory) || ...
                    ~isa(solver.solverFactory, 'function_handle')
                return
            end
            name = class(solver.solverFactory(CPPLINE.probeNetwork()));
            short = upper(name);
            if strncmp(short, 'SOLVER', 6)
                short = short(7:end);
            end
            switch short
                case 'MVA'
                    layerSolverType = jline.lang.constant.SolverType.MVA;
                case {'NC', 'COMOM'}
                    layerSolverType = jline.lang.constant.SolverType.NC;
                case {'FLD', 'FLUID'}
                    layerSolverType = jline.lang.constant.SolverType.FLUID;
                case 'CTMC'
                    layerSolverType = jline.lang.constant.SolverType.CTMC;
                case 'MAM'
                    layerSolverType = jline.lang.constant.SolverType.MAM;
                case 'SSA'
                    layerSolverType = jline.lang.constant.SolverType.SSA;
                case 'JMT'
                    layerSolverType = jline.lang.constant.SolverType.JMT;
                case 'QNS'
                    layerSolverType = jline.lang.constant.SolverType.QNS;
                case 'AUTO'
                    layerSolverType = jline.lang.constant.SolverType.AUTO;
                otherwise
                    line_error(mfilename, sprintf(['lang=''java'' runs the LQN layers under ' ...
                        'SolverAUTO, SolverCTMC, SolverFluid, SolverJMT, SolverMAM, SolverMVA, ' ...
                        'SolverNC, SolverQNS or SolverSSA; this SolverLN builds a ''%s'' layer ' ...
                        'solver, which the JAR has no layer factory for. Solve it with ' ...
                        'lang=''matlab'', or build the SolverLN with one of those layer ' ...
                        'factories.'], name));
            end
        end

        function jfun = reward_handle_to_tabulatedfun(rewardFn, sn)
            % REWARD_HANDLE_TO_TABULATEDFUN Convert a MATLAB reward function
            % handle to a Java TabulatedRewardFunction by pre-computing its
            % value over an enumerable superset of the aggregated state space.
            %
            % The domain per class r is: all per-station count vectors with
            % sum <= njobs(r) for closed classes, and per-station counts capped
            % by classcap(i,r) (or the default CTMC cutoff of 100 when
            % infinite) for open classes. The JAR reward analyzer evaluates the
            % function only on reachable stateSpaceAggr rows, which are a
            % subset of this domain; unseen states raise a descriptive error.
            %
            % @param rewardFn MATLAB reward handle @(state) or @(state, sn)
            % @param sn Network struct
            % @return jfun Java jline.lang.reward.TabulatedRewardFunction

            M = sn.nstations;
            K = sn.nclasses;

            % Per-class families of feasible per-station count vectors, with a
            % guard against combinatorial explosion (as in the PAS precompute)
            classVecs = cell(1, K);
            total = 1;
            for r = 1:K
                if isfinite(sn.njobs(r)) && sn.njobs(r) > 0
                    classVecs{r} = JLINE.reward_enum_capped(M, sn.njobs(r)*ones(1,M), sn.njobs(r));
                else
                    caps = zeros(1, M);
                    for i = 1:M
                        if sn.sched(i) == SchedStrategy.EXT
                            % Source never holds jobs in the aggregated CTMC state -- see _kb/04-networkstruct.md
                            caps(i) = 0;
                            continue;
                        end
                        c = sn.classcap(i,r);
                        if ~isfinite(c)
                            c = 100; % default CTMC cutoff (see solver_ctmc_reward)
                        end
                        caps(i) = c;
                    end
                    classVecs{r} = JLINE.reward_enum_capped(M, caps, Inf);
                end
                total = total * size(classVecs{r},1);
                if total > 5e6
                    line_error(mfilename, 'Reward pre-computation would enumerate more than 5e6 aggregated states; set finite class capacities (or smaller populations) for JLINE conversion.');
                end
            end

            % Index maps for RewardState, as in solver_ctmc_reward
            nodeToStationMap = configureDictionary('int32', 'int32');
            classToIndexMap = configureDictionary('int32', 'int32');
            for ind = 1:sn.nnodes
                if sn.isstation(ind)
                    nodeToStationMap(int32(ind)) = sn.nodeToStation(ind);
                end
            end
            for r = 1:K
                classToIndexMap(int32(r)) = r;
            end

            jfun = javaObject('jline.lang.reward.TabulatedRewardFunction');
            counts = zeros(1, K);
            for r = 1:K
                counts(r) = size(classVecs{r}, 1);
            end
            idx = ones(1, K);
            while true
                row = zeros(1, M*K);
                for r = 1:K
                    row(((1:M)-1)*K + r) = classVecs{r}(idx(r), :);
                end
                rewardState = RewardState(row, sn, nodeToStationMap, classToIndexMap);
                % Try the new single-argument API first, then the backward
                % compatible @(state, sn) signature (as in solver_ctmc_reward)
                try
                    val = rewardFn(rewardState);
                catch ME
                    try
                        val = rewardFn(row, sn);
                    catch
                        rethrow(ME);
                    end
                end
                jfun.addValue(JLINE.from_line_matrix(row), double(val));
                % Advance the mixed-radix odometer over classes
                r = 1;
                while r <= K
                    idx(r) = idx(r) + 1;
                    if idx(r) <= counts(r)
                        break;
                    end
                    idx(r) = 1;
                    r = r + 1;
                end
                if r > K
                    break;
                end
            end
        end

        function V = reward_enum_capped(M, caps, budget)
            % REWARD_ENUM_CAPPED All integer row vectors v (1 x M) with
            % 0 <= v(i) <= caps(i) and sum(v) <= budget.
            if M == 1
                hi = min(caps(1), budget);
                V = (0:hi)';
                return
            end
            V = zeros(0, M);
            hi = min(caps(1), budget);
            for n = 0:hi
                Vsub = JLINE.reward_enum_capped(M-1, caps(2:end), budget - n);
                V = [V; [n*ones(size(Vsub,1),1), Vsub]]; %#ok<AGROW>
            end
        end

        function serfun = pas_handle_to_serializablefun(handle, nclasses, cap)
            % PAS_HANDLE_TO_SERIALIZABLEFUN Convert a PAS service rate function
            % mu(c) (MATLAB handle of the ordered class list) to a Java
            % SerializableFunction by pre-computing mu over every ordered prefix
            % up to the queue capacity.
            %
            % The JAR queries mu(c) with the ordered prefix as a row vector of
            % 0-based class indices (jline.lang.state.AfterEventStation), and the
            % PrecomputedRateFunction keys on the stringified vector, so values
            % are stored under the matching 0-based key.
            %
            % @param handle   MATLAB mu(c) handle taking a 1-based ordered class list
            % @param nclasses Number of classes
            % @param cap      Station capacity (max ordered-list length)
            % @return serfun  Java PrecomputedCDFunction

            % Guard against combinatorial explosion of ordered sequences
            nseq = 0; term = 1;
            for k = 1:cap
                term = term * nclasses;
                nseq = nseq + term;
            end
            if nseq > 5e6
                line_error(mfilename, sprintf('PAS service rate pre-computation would enumerate %d ordered states (nclasses=%d, cap=%d); too large for JLINE conversion.', nseq, nclasses, cap));
            end

            serfun = jline.util.PrecomputedRateFunction(nclasses, 0.0);
            JLINE.pas_enumerate_seqs(handle, serfun, nclasses, cap, []);
        end

        function pas_enumerate_seqs(handle, serfun, nclasses, cap, prefix)
            % PAS_ENUMERATE_SEQS Recursively enumerate ordered class sequences
            % (1-based, length 1..cap) and store mu(prefix) under the matching
            % 0-based key expected by the JAR.
            if ~isempty(prefix)
                val = handle(prefix);                       % mu(c), 1-based ordered list
                keyMat = JLINE.from_line_matrix(prefix - 1); % 0-based key (1 x p)
                serfun.addValue(keyMat, double(val));
            end
            if length(prefix) >= cap
                return;
            end
            for r = 1:nclasses
                JLINE.pas_enumerate_seqs(handle, serfun, nclasses, cap, [prefix, r]);
            end
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
                % Complete state: errors NOT swallowed here (a lookup miss falls back to beta=1, i.e. unscaled)
                value = handle(currentState);
                % Convert to Java int array and add to serfun
                jstate = jline.util.matrix.Matrix(1, nclasses);
                for r = 1:nclasses
                    jstate.set(0, r-1, currentState(r));
                end
                if isscalar(value)
                    % Chain-independent beta_i(n): one scaling shared by every
                    % class, broadcast on the Java side.
                    serfun.addValue(jstate, double(value));
                else
                    % Chain-specific beta_{i,r}(n): keep the per-class vector whole (double[] overload), not the scalar one
                    serfun.addValue(jstate, double(value(:)'));
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

        function jSched = to_jline_sched_strategy(schedId)
            % Convert MATLAB SchedStrategy id to jline SchedStrategy enum
            switch schedId
                case SchedStrategy.REF
                    jSched = jline.lang.constant.SchedStrategy.REF;
                case SchedStrategy.INF
                    jSched = jline.lang.constant.SchedStrategy.INF;
                case SchedStrategy.FCFS
                    jSched = jline.lang.constant.SchedStrategy.FCFS;
                case SchedStrategy.LCFS
                    jSched = jline.lang.constant.SchedStrategy.LCFS;
                case SchedStrategy.SIRO
                    jSched = jline.lang.constant.SchedStrategy.SIRO;
                case SchedStrategy.SJF
                    jSched = jline.lang.constant.SchedStrategy.SJF;
                case SchedStrategy.LJF
                    jSched = jline.lang.constant.SchedStrategy.LJF;
                case SchedStrategy.PS
                    jSched = jline.lang.constant.SchedStrategy.PS;
                case SchedStrategy.DPS
                    jSched = jline.lang.constant.SchedStrategy.DPS;
                case SchedStrategy.GPS
                    jSched = jline.lang.constant.SchedStrategy.GPS;
                case SchedStrategy.SEPT
                    jSched = jline.lang.constant.SchedStrategy.SEPT;
                case SchedStrategy.LEPT
                    jSched = jline.lang.constant.SchedStrategy.LEPT;
                case SchedStrategy.HOL
                    % preserve HOL (exact M/G/1 priority in the JAR) rather than
                    % collapsing to FCFSPRIO (egflin approximation)
                    jSched = jline.lang.constant.SchedStrategy.HOL;
                case SchedStrategy.FCFSPRIO
                    jSched = jline.lang.constant.SchedStrategy.FCFSPRIO;
                case SchedStrategy.FORK
                    jSched = jline.lang.constant.SchedStrategy.FORK;
                case SchedStrategy.EXT
                    jSched = jline.lang.constant.SchedStrategy.EXT;
                case SchedStrategy.LCFSPR
                    jSched = jline.lang.constant.SchedStrategy.LCFSPR;
                case SchedStrategy.PSPRIO
                    jSched = jline.lang.constant.SchedStrategy.PSPRIO;
                case SchedStrategy.DPSPRIO
                    jSched = jline.lang.constant.SchedStrategy.DPSPRIO;
                case SchedStrategy.GPSPRIO
                    jSched = jline.lang.constant.SchedStrategy.GPSPRIO;
                otherwise
                    jSched = jline.lang.constant.SchedStrategy.FCFS;
            end
        end

        function tf = is_custom_handle(funCell, e, h, defaultStr)
            % IS_CUSTOM_HANDLE True if funCell{e,h} is a function handle that
            % differs from the given identity default (whitespace-insensitive
            % func2str comparison).
            tf = false;
            if isempty(funCell) || size(funCell,1) < e || size(funCell,2) < h
                return;
            end
            fh = funCell{e,h};
            if isempty(fh) || ~isa(fh, 'function_handle')
                return;
            end
            tf = ~strcmp(regexprep(func2str(fh), '\s+', ''), defaultStr);
        end

        function jwf = from_line_workflow(line_wf)
            % Convert a MATLAB Workflow to a JAR Workflow.
            jwf = javaObject('jline.lang.workflow.Workflow', java.lang.String(line_wf.getName()));
            acts = line_wf.activities;
            for a = 1:length(acts)
                act = acts{a};
                actName = java.lang.String(act.name);
                if ~isempty(act.hostDemand) && isa(act.hostDemand, 'Distribution')
                    jdist = JLINE.from_line_distribution(act.hostDemand);
                    jwf.addActivity(actName, jdist);
                else
                    jwf.addActivity(actName, 1.0);
                end
            end
            precs = line_wf.precedences;
            for p = 1:length(precs)
                prec = precs(p);
                preActs = java.util.ArrayList();
                for k = 1:length(prec.preActs)
                    preActs.add(java.lang.String(prec.preActs{k}));
                end
                postActs = java.util.ArrayList();
                for k = 1:length(prec.postActs)
                    postActs.add(java.lang.String(prec.postActs{k}));
                end
                preTypeStr = java.lang.String(ActivityPrecedenceType.toText(prec.preType));
                postTypeStr = java.lang.String(ActivityPrecedenceType.toText(prec.postType));
                if ~isempty(prec.preParams)
                    preParamsMat = JLINE.from_line_matrix(prec.preParams(:)');
                else
                    preParamsMat = javaObject('jline.util.matrix.Matrix', 0, 0);
                end
                if ~isempty(prec.postParams)
                    postParamsMat = JLINE.from_line_matrix(prec.postParams(:)');
                else
                    postParamsMat = javaObject('jline.util.matrix.Matrix', 0, 0);
                end
                jprec = javaObject('jline.lang.layered.ActivityPrecedence', ...
                    preActs, postActs, preTypeStr, postTypeStr, preParamsMat, postParamsMat);
                jwf.addPrecedence(jprec);
            end
        end

        function jenv = from_line_environment(line_env)
            % Convert a MATLAB Environment to a JAR Environment.
            E = height(line_env.envGraph.Nodes);
            jenv = javaObject('jline.lang.Environment', java.lang.String(line_env.getName()), int32(E));
            for e = 1:E
                stageName = char(line_env.envGraph.Nodes.Name{e});
                stageType = char(line_env.envGraph.Nodes.Type{e});
                stageModel = line_env.ensemble{e};
                jmodel = JLINE.from_line_network(stageModel);
                jenv.addStage(int32(e-1), java.lang.String(stageName), java.lang.String(stageType), jmodel);
            end
            if ~isempty(line_env.env)
                [Erows, Ecols] = size(line_env.env);
                for e = 1:Erows
                    for h = 1:Ecols
                        d = line_env.env{e,h};
                        if isempty(d) || isa(d, 'Disabled')
                            continue;
                        end
                        % Custom reset handles cannot marshal to JAR functional interfaces; error rather than drop silently
                        if JLINE.is_custom_handle(line_env.resetFun, e, h, '@(q)q') ...
                                || JLINE.is_custom_handle(line_env.resetEnvRatesFun, e, h, '@(originalDist,QExit,UExit,TExit)originalDist') ...
                                || JLINE.is_custom_handle(line_env.resetStateFun, e, h, '@(pi)pi')
                            line_error(mfilename, sprintf('JLINE conversion cannot marshal the custom reset function on Environment transition %d->%d (resetFun/resetEnvRatesFun/resetStateFun); use the MATLAB-native SolverENV for this model.', e, h));
                        end
                        jdist = JLINE.from_line_distribution(d);
                        jenv.addTransition(int32(e-1), int32(h-1), jdist);
                    end
                end
            end
            jenv.init();
        end

    end
end
