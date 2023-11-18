classdef JLINE

    % Example:
    % java_network = JLINE.from_line_network(network);
    % jssa = JLINE.get_solver(java_network, 'ssa');

    properties (Constant)
        jar_loc = which('jline.jar');
    end

    methods(Static)

        function model = from_line_layered_network(line_layered_network)
            sn = line_layered_network.getStruct;

            %% initialization
            model = jline.lang.LayeredNetwork(line_layered_network.getName);

            %% host processors
            P = cell(1,sn.nhosts);
            for h=1:sn.nhosts
                P{h} = jline.lang.Processor(model, sn.names{h}, sn.mult(h), jline.lang.constant.SchedStrategy.fromLINEString(sn.sched{h}));
                if sn.repl(h)~=1
                    P{h}.setReplication(sn.repl(h));
                end
            end

            %% tasks
            T = cell(1,sn.ntasks);
            for t=1:sn.ntasks
                tidx = sn.tshift+t;
                T{t} = jline.lang.Task(model, sn.names{tidx}, sn.mult(tidx), jline.lang.constant.SchedStrategy.fromLINEString(sn.sched{tidx}));
                T{t}.on(P{sn.parent(tidx)})
                if sn.repl(tidx)~=1
                    T{t}.setReplication(sn.repl(tidx));
                end
                if ~isempty(sn.think{tidx})
                    T{t}.setThinkTime(JLINE.from_line_distribution(sn.think{tidx}));
                end
            end
            %% entries
            E = cell(1,sn.nentries);
            for e=1:sn.nentries
                eidx = sn.eshift+e;
                E{e} = jline.lang.Entry(model, sn.names{eidx});
                E{e}.on(T{sn.parent(eidx)-sn.tshift});
            end

            %% activities
            A = cell(1,sn.nacts);
            for a=1:sn.nacts
                aidx = sn.ashift+a;
                tidx = sn.parent(aidx);
                onTask = tidx-sn.tshift;
                A{a} = jline.lang.Activity(model, sn.names{aidx}, JLINE.from_line_distribution(sn.hostdem{aidx}));
                A{a}.on(T{onTask});

                boundTo = find(sn.graph((sn.eshift+1):(sn.eshift+sn.nentries),aidx));

                if ~isempty(boundTo)
                    A{a}.boundTo(E{boundTo});
                end

                if sn.sched{tidx} ~= SchedStrategy.ID_REF % ref tasks don't reply
                    repliesTo = find(sn.replygraph(aidx,:))-sn.eshift; % index of entry
                    if ~isempty(repliesTo)
                        if ~sn.isref(sn.parent(sn.eshift+repliesTo))
                            A{a}.repliesTo(E{repliesTo});
                        end
                    end
                end

                if ~isempty(sn.callpair)
                    cidxs = find(sn.callpair(:,1)==aidx);
                    calls = sn.callpair(:,2);
                    for c=cidxs(:)'
                        switch sn.calltype(c)
                            case CallType.ID_SYNC
                                A{a}.synchCall(E{calls(c)-sn.eshift},sn.callproc{c}.getMean);
                            case CallType.ID_ASYNC
                                A{a}.asynchCall(E{calls(c)-sn.eshift},sn.callproc{c}.getMean);
                        end
                    end
                end

            end

            %% think times
            for h=1:sn.nhosts
                if ~isempty(sn.think{h})
                    switch sn.think{h}.name
                        case 'Immediate'
                            P{h}.setThinkTime(jline.lang.distributions.Immediate);
                        case 'Exp'
                            P{h}.setThinkTime(jline.lang.distributions.Exp(1/sn.hostdem{aidx}.getMean));
                        case 'Erlang'
                            P{h}.setThinkTime(jline.lang.distributions.Erlang.fitMeanAndSCV(sn.hostdem{aidx}.getMean,sn.hostdem{aidx}.getSCV));
                        case 'HyperExp'
                            P{h}.setThinkTime(jline.lang.distributions.HyperExp(sn.hostdem{aidx}.getParam(1).paramValue, sn.hostdem{aidx}.getParam(2).paramValue, sn.hostdem{aidx}.getParam(3).paramValue));
                        case 'Coxian'
                            P{h}.setThinkTime(jline.lang.distributions.Coxian(sn.hostdem{aidx}.getParam(1).paramValue, sn.hostdem{aidx}.getParam(2).paramValue));
                        case 'APH'
                            P{h}.setThinkTime(jline.lang.distributions.APH(sn.hostdem{aidx}.getParam(1).paramValue, sn.hostdem{aidx}.getParam(2).paramValue));
                        otherwise
                            line_error(mfilename,sprintf('JLINE conversion does not support the %d distribution yet.',sn.hostdem{aidx}.name));
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
                        if full(sn.actpretype(aidx)) == ActivityPrecedenceType.ID_PRE_SEQ && full(sn.actposttype(bidx)) == ActivityPrecedenceType.ID_POST_SEQ
                            T{tidx-sn.tshift}.addPrecedence(jline.lang.ActivityPrecedence.Sequence(sn.names{aidx}, sn.names{bidx}));
                        end
                    end
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%% translated up to here

            %% Loop precedences (POST_LOOP)
            precMarker = 0;
            for ai = 1:sn.nacts
                aidx = sn.ashift + ai;
                tidx = sn.parent(aidx);
                % for all successors
                for bidx=find(sn.graph(aidx,:))
                    if bidx > sn.ashift % ignore precedence between entries and activities
                        % Loop pattern (POST_LOOP)
                        if full(sn.actposttype(bidx)) == ActivityPrecedenceType.ID_POST_LOOP
                            if precMarker == 0 % start a new loop
                                precMarker = aidx-sn.ashift;
                                precActs = '';
                            else
                                if isempty(precActs)
                                    precActs = java.util.ArrayList();
                                    precActs.add(sprintf("%s", sn.names{bidx}));
                                else
                                    precActs.add(sprintf("%s", sn.names{bidx}));
                                end
                                if aidx ~= bidx % loop end reached
                                    counts = 1/sn.graph(aidx,bidx);
                                    T{tidx-sn.tshift}.addPrecedence(jline.lang.ActivityPrecedence.Loop(sn.names{precMarker+sn.ashift}, precActs, jline.util.Matrix(counts)));
                                    precMarker = 0;
                                end
                            end
                        end
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
                        if full(sn.actposttype(bidx)) == ActivityPrecedenceType.ID_POST_OR
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
                    probsMatrix = jline.util.Matrix(1,length(probs));
                    for i=1:length(probs)
                        probsMatrix.set(0,i-1,probs(i));
                    end
                    T{tidx-sn.tshift}.addPrecedence(jline.lang.ActivityPrecedence.OrFork(sn.names{precMarker+sn.ashift}, precActs, probsMatrix));
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
                        if full(sn.actposttype(bidx)) == ActivityPrecedenceType.ID_POST_AND
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
                    T{tidx-sn.tshift}.addPrecedence(jline.lang.ActivityPrecedence.AndFork(sn.names{precMarker+sn.ashift}, postActs));
                    precMarker = 0;
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
                        if full(sn.actpretype(aidx)) == ActivityPrecedenceType.ID_PRE_OR
                            if precMarker == 0 % start a new orjoin
                                precActs = java.util.ArrayList();
                                precMarker = bidx-sn.ashift;
                                precActs.add(sprintf("%s", sn.names{bidx}));
                            else
                                precActs.add(sprintf("%s", sn.names{bidx}));
                            end
                        end
                    end
                end
                if precMarker > 0
                    T{tidx-sn.tshift}.addPrecedence(jline.lang.ActivityPrecedence.OrJoin(precActs, sn.names{precMarker+sn.ashift}));
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
                        if full(sn.actpretype(aidx)) == ActivityPrecedenceType.ID_PRE_AND
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
                    T{tidx-sn.tshift}.addPrecedence(jline.lang.ActivityPrecedence.AndJoin(precActs, sn.names{precMarker+sn.ashift}, length(precActs)));
                    precMarker = 0;
                end
            end

        end

        function java_dist = from_line_distribution(line_dist)
            if isa(line_dist, 'Exp')
                java_dist = jline.lang.distributions.Exp(line_dist.getParam(1).paramValue);
            elseif isa(line_dist, "APH")
                alpha = line_dist.getParam(1).paramValue;
                T = line_dist.getParam(2).paramValue;
                jline_alpha = java.util.ArrayList();
                for i = 1:length(alpha)
                    jline_alpha.add(alpha(i));
                end
                jline_T = JLINE.matrix_to_jlinematrix(T);
                java_dist = jline.lang.distributions.APH(jline_alpha, jline_T);
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
                java_dist = jline.lang.distributions.Coxian(jline_mu, jline_phi);
            elseif isa(line_dist, 'Det')
                java_dist = jline.lang.distributions.Det(line_dist.getParam(1).paramValue);
            elseif isa(line_dist, 'Erlang')
                java_dist = jline.lang.distributions.Erlang(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Gamma')
                java_dist = jline.lang.distributions.Gamma(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, "HyperExp")
                java_dist = jline.lang.distributions.HyperExp(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue, line_dist.getParam(3).paramValue);
            elseif isa(line_dist, 'LogNormal')
                java_dist = jline.lang.distributions.LogNormal(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Pareto')
                java_dist = jline.lang.distributions.Pareto(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, "PH")
                alpha = line_dist.getParam(1).paramValue;
                T = line_dist.getParam(2).paramValue;
                jline_alpha = java.util.ArrayList();
                for i = 1:length(alpha)
                    jline_alpha.add(alpha(i));
                end
                jline_T = JLINE.matrix_to_jlinematrix(T);
                java_dist = jline.lang.distributions.PH(jline_alpha, jline_T);
            elseif isa(line_dist, 'Uniform')
                java_dist = jline.lang.distributions.Uniform(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Weibull')
                java_dist = jline.lang.distributions.Weibull(line_dist.getParam(1).paramValue, line_dist.getParam(2).paramValue);
            elseif isa(line_dist, 'Immediate')
                java_dist = jline.lang.distributions.Immediate();
            elseif isempty(line_dist) || isa(line_dist, 'Disabled')
                java_dist = jline.lang.distributions.Disabled();
                return;
            elseif isa(line_dist, 'Replayer')
                line_error(mfilename,'The Replayer distribution is not yet supported by JLINE.');
            else
                line_error(mfilename,'Distribution not supported by JLINE.');
            end
        end

        function matlab_dist = from_jline_distribution(java_dist)
            if isa(java_dist, 'jline.lang.distributions.Exp')
                matlab_dist = Exp(java_dist.getRate());
            elseif isa(java_dist, 'jline.lang.distributions.Det')
                matlab_dist = Det(java_dist.getParam(1).getValue);
            elseif isa(java_dist, 'jline.lang.distributions.Erlang')
                matlab_dist = Erlang(java_dist.getRate(),java_dist.getNumberOfPhases());
            elseif isa(java_dist, 'jline.lang.distributions.Gamma')
                matlab_dist = Gamma(java_dist.getParam(1).getValue,java_dist.getParam().getValue);
            elseif isa(java_dist, 'jline.lang.distributions.HyperExp')
                matlab_dist = HyperExp(java_dist.getParam(1).getValue, java_dist.getParam(2).getValue, java_dist.getParam(3).getValue);
            elseif isa(java_dist, 'jline.lang.distributions.LogNormal')
                matlab_dist = LogNormal(java_dist.getParam(1).getValue,java_dist.getParam().getValue);
            elseif isa(java_dist, 'jline.lang.distributions.Pareto')
                matlab_dist = Pareto(java_dist.getParam(1).getValue,java_dist.getParam().getValue);
            elseif isa(java_dist, 'jline.lang.distributions.Uniform')
                matlab_dist = Uniform(java_dist.getParam(1).getValue,java_dist.getParam().getValue);
            elseif isa(java_dist, 'jline.lang.distributions.Weibull')
                matlab_dist = Weibull(java_dist.getParam(1).getValue,java_dist.getParam().getValue);
            elseif isa(java_dist, 'jline.lang.distributions.APH')
                java_alpha = java_dist.getParam(1).getValue;
                alpha = zeros(1, java_alpha.size);
                for i = 1:java_alpha.size
                    alpha(i) = java_alpha.get(i-1);
                end
                matlab_dist = APH(alpha, JLINE.jlinematrix_to_matrix(java_dist.getParam(2).getValue));
            elseif isa(java_dist, 'jline.lang.distributions.PH')
                java_alpha = java_dist.getParam(1).getValue;
                alpha = zeros(1, java_alpha.size);
                for i = 1:java_alpha.size
                    alpha(i) = java_alpha.get(i-1);
                end
                matlab_dist = PH(alpha, JLINE.jlinematrix_to_matrix(java_dist.getParam(2).getValue));
            elseif isa(java_dist, 'jline.lang.distributions.Coxian')
                if java_dist.getNumberOfPhases == 2
                    matlab_dist = Coxian(java_dist.getParam(1).getValue.get(0), java_dist.getParam(1).getValue.get(1), java_dist.getParam(2).getValue.get(0));
                else
                    java_mu = java_dist.getParam(1).getValue;
                    java_phi = java_dist.getParam(2).getValue;
                    mu = zeros(1, java_mu.size);
                    phi = zeros(1, java_phi.size);
                    for i = 1:java_mu.size
                        mu(i) = java_mu.get(i-1);
                    end
                    for i = 1:java_phi.size
                        phi(i) = java_phi.get(i-1);
                    end
                    matlab_dist = Coxian(mu, phi);
                end
            elseif isa(java_dist, 'jline.lang.distributions.Immediate')
                matlab_dist = Immediate();
            elseif isa(java_dist, 'jline.lang.distributions.Disabled')
                return;
            else
                line_error(mfilename,'Distribution not supported by JLINE.');
            end
        end

        function set_csMatrix(line_node, java_node)
            nClasses = length(line_node.server.classes);
            csMatrix = jline.util.Matrix(nClasses, nClasses);
            for i = 1:nClasses
                for j = 1:nClasses
                    csMatrix.set(i-1, j-1, line_node.server.csFun(i,j,0,0));
                end
            end
            java_node.setClassSwitchingMatrix(csMatrix);
        end

        function set_service(line_node, java_node, job_classes)
            if (isa(line_node, 'Sink') || isa(line_node, 'Router') || isa(line_node, 'ClassSwitch') || isa(line_node, 'Fork') || isa(line_node, 'Join'))
                return;
            end

            for n = 1 : length(job_classes)
                if (isa(line_node, 'Queue') || isa(line_node, 'Delay'))
                    matlab_dist = line_node.getService(job_classes{n});
                elseif (isa(line_node, 'Source'))
                    matlab_dist = line_node.getArrivalProcess(n);
                else
                    line_error(mfilename,'Node not supported by JLINE.');
                end
                service_dist = JLINE.from_line_distribution(matlab_dist);

                if (isa(line_node,'Queue') || isa(line_node, 'Delay'))
                    java_node.setService(java_node.getModel().getClasses().get(n-1), service_dist, line_node.schedStrategyPar(n));
                elseif (isa(line_node, 'Source'))
                    java_node.setArrival(java_node.getModel().getClasses().get(n-1), service_dist);
                end
            end
        end

        function set_line_service(jline_node, line_node, job_classes, line_classes)
            if (isa(line_node,'Sink')) || isa(line_node, 'ClassSwitch')
                return;
            end
            for n = 1:job_classes.size()
                java_dist = jline_node.getServiceProcess(job_classes.get(n-1));
                matlab_dist = JLINE.from_jline_distribution(java_dist);

                if (isa(line_node, 'Queue') || isa(line_node, 'Delay'))
                    line_node.setService(line_classes{n}, matlab_dist);
                elseif (isa(line_node, 'Source'))
                    line_node.setArrival(line_classes{n}, matlab_dist);
                else
                    line_error(mfilename,'Node not supported by JLINE.');
                end
            end
        end

        function node_object = from_line_node(line_node, java_network, ~, forkNode)
            if isa(line_node, 'Delay')
                node_object = jline.lang.nodes.Delay(java_network, line_node.getName);
            elseif isa(line_node, 'Queue')
                switch SchedStrategy.toId(line_node.schedStrategy)
                    case SchedStrategy.ID_INF
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.INF);
                    case SchedStrategy.ID_FCFS
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.FCFS);
                    case SchedStrategy.ID_LCFS
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.LCFS);
                    case SchedStrategy.ID_SIRO
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.SIRO);
                    case SchedStrategy.ID_SJF
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.SJF);
                    case SchedStrategy.ID_LJF
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.LJF);
                    case SchedStrategy.ID_PS
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.PS);
                    case SchedStrategy.ID_DPS
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.DPS);
                    case SchedStrategy.ID_GPS
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.GPS);
                    case SchedStrategy.ID_SEPT
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.SEPT);
                    case SchedStrategy.ID_LEPT
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.LEPT);
                    case SchedStrategy.ID_HOL
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.HOL);
                    case SchedStrategy.ID_FORK
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.FORK);
                    case SchedStrategy.ID_EXT
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.EXT);
                    case SchedStrategy.ID_REF
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.REF);
                    case SchedStrategy.ID_LCFSPR
                        node_object = jline.lang.nodes.Queue(java_network, line_node.getName, jline.lang.constant.SchedStrategy.LCFSPR);
                end
                if isinf(line_node.getNumberOfServers)
                    node_object.setNumberOfServers(java.lang.Integer.MAX_VALUE);
                else
                    node_object.setNumberOfServers(line_node.getNumberOfServers);
                end
                if ~isempty(line_node.lldScaling)
                    node_object.setLoadDependence(JLINE.matrix_to_jlinematrix(line_node.lldScaling));
                end
            elseif isa(line_node, 'Source')
                node_object = jline.lang.nodes.Source(java_network, line_node.getName);
            elseif isa(line_node, 'Sink')
                node_object = jline.lang.nodes.Sink(java_network, line_node.getName);
            elseif isa(line_node, 'Router')
                node_object = jline.lang.nodes.Router(java_network, line_node.getName);
            elseif isa(line_node, 'ClassSwitch')
                node_object = jline.lang.nodes.ClassSwitch(java_network, line_node.getName);
            elseif isa(line_node, 'Fork')
                node_object = jline.lang.nodes.Fork(java_network);
            elseif isa(line_node, 'Join')
                node_object = jline.lang.nodes.Join(java_network, line_node.name, forkNode);
            else
                line_error(mfilename,'Node not supported by JLINE.');
            end
        end

        function node_object = from_jline_node(jline_node, line_network, job_classes)
            if isa(jline_node, 'jline.lang.nodes.Delay')
                node_object = Delay(line_network, jline_node.getName.toCharArray');
            elseif isa(jline_node, 'jline.lang.nodes.Queue')
                schedStrategy = jline_node.getSchedStrategy;
                switch schedStrategy.name().toCharArray'
                    case 'INF'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.INF);
                    case 'FCFS'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.FCFS);
                    case 'LCFS'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.LCFS);
                    case 'SIRO'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.SIRO);
                    case 'SJF'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.SJF);
                    case 'LJF'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.LJF);
                    case 'PS'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.PS);
                    case 'DPS'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.DPS);
                    case 'GPS'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.GPS);
                    case 'SEPT'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.SEPT);
                    case 'LEPT'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.LEPT);
                    case 'HOL'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.HOL);
                    case 'FORK'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.FORK);
                    case 'EXT'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.EXT);
                    case 'REF'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.REF);
                    case 'LCFSPR'
                        node_object = Queue(line_network, jline_node.getName.toCharArray', SchedStrategy.LCFSPR);
                end
                node_object.setNumberOfServers(jline_node.getNumberOfServers);
                if ~isempty(JLINE.jlinematrix_to_matrix(jline_node.getLimitedLoadDependence))
                    node_object.setLoadDependence(JLINE.jlinematrix_to_matrix(jline_node.getLimitedLoadDependence));
                end
            elseif isa(jline_node, 'jline.lang.nodes.Source')
                node_object = Source(line_network, jline_node.getName.toCharArray');
            elseif isa(jline_node, 'jline.lang.nodes.Sink')
                node_object = Sink(line_network, jline_node.getName.toCharArray');
            elseif isa(jline_node, 'jline.lang.nodes.Router')
                node_object = Router(line_network, jline_node.getName.toCharArray');
            elseif isa(jline_node, 'jline.lang.nodes.ClassSwitch')
                nClasses = job_classes.size;
                csMatrix = zeros(nClasses, nClasses);
                for r = 1:nClasses
                    for s = 1:nClasses
                        csMatrix(r,s) = jline_node.getServer.applyCsFun(r-1,s-1);
                    end
                end
                node_object = ClassSwitch(line_network, jline_node.getName.toCharArray', csMatrix);
            else
                line_error(mfilename,'Node not supported by JLINE.');
            end
        end

        function node_class = from_line_class(line_class, java_network)
            if isa(line_class, 'OpenClass')
                node_class = jline.lang.OpenClass(java_network, line_class.getName, line_class.priority);
            elseif isa(line_class, 'SelfLoopingClass')
                node_class = jline.lang.SelfLoopingClass(java_network, line_class.getName, line_class.population, java_network.getNodeByName(line_class.refstat.getName), line_class.priority);
            elseif isa(line_class, 'ClosedClass')
                node_class = jline.lang.ClosedClass(java_network, line_class.getName, line_class.population, java_network.getNodeByName(line_class.refstat.getName), line_class.priority);
            else
                line_error(mfilename,'Class type not supported by JLINE.');
            end
        end

        function node_class = from_jline_class(java_class, line_network)
            if isa(java_class, 'jline.lang.OpenClass')
                node_class = OpenClass(line_network, java_class.getName.toCharArray', java_class.getPriority);
            elseif isa(java_class, 'jline.lang.SelfLoopingClass')
                node_class = SelfLoopingClass(line_network, java_class.getName.toCharArray', java_class.getNumberOfJobs, line_network.getNodeByName(java_class.getRefstat.getName), java_class.getPriority);
            elseif isa(java_class, 'jline.lang.ClosedClass')
                node_class = ClosedClass(line_network, java_class.getName.toCharArray', java_class.getNumberOfJobs, line_network.getNodeByName(java_class.getRefstat.getName), java_class.getPriority);
            else
                line_error(mfilename,'Class type not supported by JLINE.');
            end
        end

        function from_line_links(line_network, network_object)
            connections = line_network.getConnectionMatrix();
            [m, ~] = size(connections);
            nodes = network_object.getNodes();
            classes = network_object.getClasses();
            n_classes = classes.size();
            routing_matrix = jline.lang.RoutingMatrix(network_object, classes, nodes);
            line_nodes = line_network.getNodes;
            % [ ] Update to consider different weights/routing for classes
            for i = 1:m
                line_node = line_nodes{i};
                for k = 1:n_classes                    
                    output_strat = line_node.output.outputStrategy{k};
                    switch output_strat{2}
                        case 'Disabled'
                            nodes.get(i-1).setRouting(classes.get(k-1),jline.lang.constant.RoutingStrategy.DISABLED);
                        case 'Random'
                            nodes.get(i-1).setRouting(classes.get(k-1),jline.lang.constant.RoutingStrategy.RAND);
                        case 'Probabilities'
                            if length(output_strat) >= 3
                                probabilities = output_strat{3};
                                for j = 1:length(probabilities)
                                    dest_idx = probabilities{j}{1}.index;
                                    if (connections(i, dest_idx) ~= 0)
                                        routing_matrix.addConnection(nodes.get(i-1), nodes.get(dest_idx-1), classes.get(k-1), probabilities{j}{2});
                                    end
                                end
                            end
                        otherwise
                            %keyboard
                            line_warning(mfilename, sprintf('''%s'' routing strategy not supported by JLINE, setting as Disabled.\n',output_strat{2}));
                            nodes.get(i-1).setRouting(classes.get(k-1),jline.lang.constant.RoutingStrategy.DISABLED);
                    end
                end
            end
            network_object.link(routing_matrix);

            %Align the sn.rtorig be the same
            sn = network_object.getStruct(true);
            rtorig = java.util.HashMap();
            if ~isempty(line_network.sn.rtorig)
                if iscell(line_network.sn.rtorig)
                    for r = 1:n_classes
                        sub_rtorig = java.util.HashMap();
                        for s = 1:n_classes
                            sub_rtorig.put(classes.get(s-1), JLINE.matrix_to_jlinematrix(line_network.sn.rtorig{r,s}));
                        end
                        rtorig.put(classes.get(r-1), sub_rtorig);
                    end
                end
            end
            sn.rtorig = rtorig;
        end

        function line_network = from_jline_links(line_network, java_network)
            P = line_network.initRoutingMatrix;
            java_nodes = java_network.getNodes();
            java_classes = java_network.getClasses();
            n_classes = java_classes.size();
            n_nodes = java_nodes.size();

            for n = 1 : n_nodes
                java_node = java_nodes.get(n-1);
                output_strategies = java_node.getOutputStrategies();
                n_strategies = output_strategies.size();
                for m = 1 : n_strategies
                    output_strat = output_strategies.get(m-1);
                    dest = output_strat.getDestination();
                    in_idx = java_network.getNodeIndex(java_node)+1;
                    out_idx = java_network.getNodeIndex(dest)+1;
                    if n_classes == 1
                        P{1}(in_idx,out_idx) = output_strat.getProbability();
                    else
                        strat_class = output_strat.getJobClass();
                        class_idx = java_network.getJobClassIndex(strat_class)+1;
                        P{class_idx,class_idx}(in_idx,out_idx) = output_strat.getProbability();
                    end
                end
            end
            line_network.link(P);

            %Align the sn.rtorig be the same (Assume Java network is
            %created by calling Network.link)
            sn = java_network.getStruct;
            rtorig = cell(n_classes, n_classes);
            for r = 1:n_classes
                for s = 1:n_classes
                    rtorig{r,s} = JLINE.jlinematrix_to_matrix(sn.rtorig.get(java_classes.get(r-1)).get(java_classes.get(s-1)));
                end
            end
            line_network.sn.rtorig = rtorig;
        end

        function [network_object] = from_line_network(line_network)
            w = warning;
            warning('off');

            warning(w);
            routing_probs = line_network.getRoutingStrategies;
            for n = 1 : length(routing_probs)
                if (routing_probs(n) ~= RoutingStrategy.ID_RAND) && (routing_probs(n) ~= RoutingStrategy.ID_PROB)
                    line_error(mfilename,'Routing strategy not supported by JLINE integration script.');
                end
            end
            network_object = jline.lang.Network(line_network.getName);
            network_nodes = line_network.getNodes;
            job_classes = line_network.classes;

            java_nodes = cell(1,length(network_nodes));
            java_classes = cell(1,length(job_classes));

            for n = 1 : length(network_nodes)
                if isa(network_nodes{n}, 'Join')
                    java_nodes{n} = JLINE.from_line_node(network_nodes{n}, network_object, job_classes, java_nodes{network_nodes{n}.joinOf.index});
                else
                    java_nodes{n} = JLINE.from_line_node(network_nodes{n}, network_object, job_classes);
                end
            end

            for n = 1 : length(job_classes)
                java_classes{n} = JLINE.from_line_class(job_classes{n}, network_object);
            end

            for n = 1: length(java_nodes)
                JLINE.set_service(network_nodes{n}, java_nodes{n}, job_classes);
            end

            for n = 1: length(java_nodes)
                if isa(network_nodes{n},"ClassSwitch")
                    JLINE.set_csMatrix(network_nodes{n}, java_nodes{n});
                elseif isa(network_nodes{n},"Join")
                    java_nodes{n}.initJoinJobClasses();
                end
            end

            % Assume JLINE and LINE network are both created via link
            JLINE.from_line_links(line_network, network_object);
        end

        function line_network = jline_to_line(java_network)
            %javaaddpath(jar_loc);
            line_network = Network(java_network.getName);
            network_nodes = java_network.getNodes;
            job_classes = java_network.getClasses;

            line_nodes = cell(network_nodes.size);
            line_classes = cell(job_classes.size);


            for n = 1 : network_nodes.size
                if ~isa(network_nodes.get(n-1), 'jline.lang.nodes.ClassSwitch')
                    line_nodes{n} = JLINE.from_jline_node(network_nodes.get(n-1), line_network, job_classes);
                end
            end

            for n = 1 : job_classes.size
                line_classes{n} = JLINE.from_jline_class(job_classes.get(n-1), line_network);
            end

            for n = 1 : network_nodes.size
                if isa(network_nodes.get(n-1), 'jline.lang.nodes.ClassSwitch')
                    line_nodes{n} = JLINE.from_jline_node(network_nodes.get(n-1), line_network, job_classes);
                end
            end

            for n = 1 : network_nodes.size
                JLINE.set_line_service(network_nodes.get(n-1), line_nodes{n}, job_classes, line_classes);
            end

            % Assume JLINE and LINE network are both created via link
            line_network = JLINE.from_jline_links(line_network, java_network);
        end

        function matrix = jlinematrix_to_matrix(jline_matrix)
            if isempty(jline_matrix)
                matrix = [];
            else
                matrix = zeros(jline_matrix.getNumRows(), jline_matrix.getNumCols());
                for row = 1:jline_matrix.getNumRows()
                    for col = 1:jline_matrix.getNumCols()
                        matrix(row, col) = jline_matrix.get(row-1, col-1);
                    end
                end
            end
        end

        function jline_matrix = matrix_to_jlinematrix(matrix)
            [rows, cols] = size(matrix);
            jline_matrix = jline.util.Matrix(rows, cols);
            for row = 1:rows
                for col = 1:cols
                    if matrix(row,col) ~= 0
                        jline_matrix.set(row-1, col-1, matrix(row, col));
                    end
                end
            end
        end

        function sn = from_jline_struct(java_network, java_sn)
            %lst, rtfun and cdscaling are not implemented
            %Due to the transformation of Java lambda to matlab function
            if nargin<2
                java_sn = java_network.getStruct(false);
            end
            java_classes = java_network.getClasses;
            java_nodes = java_network.getNodes;
            java_stations = java_network.getStations;
            sn = NetworkStruct();

            sn.nnodes = java_sn.nnodes;
            sn.nclasses = java_sn.nclasses;
            sn.nclosedjobs = java_sn.nclosedjobs;
            sn.nstations = java_sn.nstations;
            sn.nstateful = java_sn.nstateful;
            sn.nchains = java_sn.nchains;

            sn.refstat = JLINE.jlinematrix_to_matrix(java_sn.refstat) + 1;
            sn.njobs = JLINE.jlinematrix_to_matrix(java_sn.njobs);
            sn.nservers = JLINE.jlinematrix_to_matrix(java_sn.nservers);
            sn.connmatrix = JLINE.jlinematrix_to_matrix(java_sn.connmatrix);
            sn.scv = JLINE.jlinematrix_to_matrix(java_sn.scv);
            sn.isstation = logical(JLINE.jlinematrix_to_matrix(java_sn.isstation));
            sn.isstateful = logical(JLINE.jlinematrix_to_matrix(java_sn.isstateful));
            sn.isstatedep = logical(JLINE.jlinematrix_to_matrix(java_sn.isstatedep));
            sn.nodeToStateful = JLINE.jlinematrix_to_matrix(java_sn.nodeToStateful)+1;
            sn.nodeToStateful(sn.nodeToStateful==0) = nan;
            sn.nodeToStation = JLINE.jlinematrix_to_matrix(java_sn.nodeToStation)+1;
            sn.nodeToStation(sn.nodeToStation==0) = nan;
            sn.stationToNode = JLINE.jlinematrix_to_matrix(java_sn.stationToNode)+1;
            sn.stationToNode(sn.stationToNode==0) = nan;
            sn.stationToStateful = JLINE.jlinematrix_to_matrix(java_sn.stationToStateful)+1;
            sn.stationToStateful(sn.stationToStateful==0) = nan;
            sn.statefulToNode = JLINE.jlinematrix_to_matrix(java_sn.statefulToNode)+1;
            sn.statefulToNode(sn.statefulToNode==0) = nan;
            sn.rates = JLINE.jlinematrix_to_matrix(java_sn.rates);
            sn.classprio = JLINE.jlinematrix_to_matrix(java_sn.classprio);
            sn.phases = JLINE.jlinematrix_to_matrix(java_sn.phases);
            sn.phasessz = JLINE.jlinematrix_to_matrix(java_sn.phasessz);
            sn.phaseshift = JLINE.jlinematrix_to_matrix(java_sn.phaseshift);
            sn.schedparam = JLINE.jlinematrix_to_matrix(java_sn.schedparam);
            sn.chains = logical(JLINE.jlinematrix_to_matrix(java_sn.chains));
            sn.rt = JLINE.jlinematrix_to_matrix(java_sn.rt);
            sn.nvars = JLINE.jlinematrix_to_matrix(java_sn.nvars);
            sn.rtnodes = JLINE.jlinematrix_to_matrix(java_sn.rtnodes);
            sn.csmask = logical(JLINE.jlinematrix_to_matrix(java_sn.csmask));
            sn.isslc = logical(JLINE.jlinematrix_to_matrix(java_sn.isslc));
            sn.cap = JLINE.jlinematrix_to_matrix(java_sn.cap);
            sn.classcap = JLINE.jlinematrix_to_matrix(java_sn.classcap);
            sn.refclass = JLINE.jlinematrix_to_matrix(java_sn.refclass)+1;
            sn.lldscaling = JLINE.jlinematrix_to_matrix(java_sn.lldscaling);

            if ~isempty(java_sn.cdscaling)
                %Not implemented since related to lambda function
            else
                sn.cdscaling = cell(sn.nstations, 0);
            end

            if ~isempty(java_sn.nodetypes)
                sn.nodetype = zeros(sn.nnodes, 1);
                for i = 1:java_sn.nodetypes.size
                    nodetype = java_sn.nodetypes.get(i-1);
                    switch nodetype.name().toCharArray'
                        case 'Queue'
                            sn.nodetype(i) = NodeType.ID_QUEUE;
                        case 'Delay'
                            sn.nodetype(i) = NodeType.ID_DELAY;
                        case 'Source'
                            sn.nodetype(i) = NodeType.ID_SOURCE;
                        case 'Sink'
                            sn.nodetype(i) = NodeType.ID_SINK;
                        case 'Join'
                            sn.nodetype(i) = NodeType.ID_JOIN;
                        case 'Fork'
                            sn.nodetype(i) = NodeType.ID_FORK;
                        case 'ClassSwitch'
                            sn.nodetype(i) = NodeType.ID_CLASSSWITCH;
                        case 'Logger'
                            sn.nodetype(i) = NodeType.ID_LOGGER;
                        case 'Cache'
                            sn.nodetype(i) = NodeType.ID_CACHE;
                        case 'Place'
                            sn.nodetype(i) = NodeType.ID_PLACE;
                        case 'Transition'
                            sn.nodetype(i) = NodeType.ID_TRANSITION;
                        case 'Router'
                            sn.nodetype(i) = NodeType.ID_ROUTER;
                    end
                end
            else
                sn.nodetype = [];
            end

            if ~isempty(java_sn.classnames)
                for i = 1:java_sn.classnames.size
                    sn.classnames(i,1) = java_sn.classnames.get(i-1);
                end
            else
                sn.classnames = [];
            end

            if ~isempty(java_sn.nodenames)
                for i = 1:java_sn.nodenames.size
                    sn.nodenames(i,1) = java_sn.nodenames.get(i-1);
                end
            else
                sn.nodenames = [];
            end

            if ~isempty(java_sn.rtorig) && java_sn.rtorig.size()>0
                sn.rtorig = cell(sn.nclasses, sn.nclasses);
                for r = 1:sn.nclasses
                    for s = 1:sn.nclasses
                        sn.rtorig{r,s} = JLINE.jlinematrix_to_matrix(java_sn.rtorig.get(java_classes.get(r-1)).get(java_classes.get(s-1)));
                    end
                end
            else
                sn.rtorig = [];
            end

            if ~isempty(java_sn.state)
                sn.state = cell(sn.nstations, 1);
                for i = 1:sn.nstations
                    sn.state{i} = JLINE.jlinematrix_to_matrix(java_sn.state.get(java_stations.get(i-1)));
                end
            else
                sn.state = {};
            end

            if ~isempty(java_sn.space)
                sn.space = cell(sn.nstations, 1);
                for i = 1:sn.nstations
                    sn.space{i} = JLINE.jlinematrix_to_matrix(java_sn.space.get(java_stations.get(i-1)));
                end
            else
                sn.space = {};
            end

            if ~isempty(java_sn.routing)
                sn.routing = zeros(sn.nnodes, sn.nclasses);
                for i = 1:sn.nnodes
                    for j = 1:sn.nclasses
                        routingStrategy = java_sn.routing.get(java_nodes.get(i-1)).get(java_classes.get(j-1));
                        switch routingStrategy.name().toCharArray'
                            case 'PROB'
                                sn.routing(i,j) = RoutingStrategy.ID_PROB;
                            case 'RAND'
                                sn.routing(i,j) = RoutingStrategy.ID_RAND;
                            case 'RROBIN'
                                sn.routing(i,j) = RoutingStrategy.ID_RROBIN;
                            case 'WRROBIN'
                                sn.routing(i,j) = RoutingStrategy.ID_WRROBIN;
                            case 'JSQ'
                                sn.routing(i,j) = RoutingStrategy.ID_JSQ;
                            case 'DISABLED'
                                sn.routing(i,j) = RoutingStrategy.ID_DISABLED;
                            case 'FIRING'
                                sn.routing(i,j) = RoutingStrategy.ID_FIRING;
                            case 'KCHOICES'
                                sn.routing(i,j) = RoutingStrategy.ID_KCHOICES;
                        end
                    end
                end
            else
                sn.routing = [];
            end

            if ~isempty(java_sn.proctype)
                sn.procid = zeros(sn.nstations, sn.nclasses);
                for i = 1:sn.nstations
                    for j = 1:sn.nclasses
                        processType = java_sn.proctype.get(java_stations.get(i-1)).get(java_classes.get(j-1));
                        %Only EXP, ERLANG, HYPEREXP, APH, COXIAN, IMMEDIATE
                        %DISABLED are implemented in JLINE.
                        switch processType.name.toCharArray'
                            case 'EXP'
                                sn.procid(i,j) = ProcessType.ID_EXP;
                            case 'ERLANG'
                                sn.procid(i,j) = ProcessType.ID_ERLANG;
                            case 'HYPEREXP'
                                sn.procid(i,j) = ProcessType.ID_HYPEREXP;
                            case 'PH'
                                sn.procid(i,j) = ProcessType.ID_PH;
                            case 'APH'
                                sn.procid(i,j) = ProcessType.ID_APH;
                            case 'MAP'
                                sn.procid(i,j) = ProcessType.ID_MAP;
                            case 'UNIFORM'
                                sn.procid(i,j) = ProcessType.ID_UNIFORM;
                            case 'DET'
                                sn.procid(i,j) = ProcessType.ID_DET;
                            case 'COXIAN'
                                sn.procid(i,j) = ProcessType.ID_COXIAN;
                            case 'GAMMA'
                                sn.procid(i,j) = ProcessType.ID_GAMMA;
                            case 'PARETO'
                                sn.procid(i,j) = ProcessType.ID_PARETO;
                            case 'WEIBULL'
                                sn.procid(i,j) = ProcessType.ID_WEIBULL;
                            case 'LOGNORMAL'
                                sn.procid(i,j) = ProcessType.ID_LOGNORMAL;
                            case 'MMPP2'
                                sn.procid(i,j) = ProcessType.ID_MMPP2;
                            case 'REPLAYER'
                                sn.procid(i,j) = ProcessType.ID_REPLAYER;
                            case 'TRACE'
                                sn.procid(i,j) = ProcessType.ID_TRACE;
                            case 'IMMEDIATE'
                                sn.procid(i,j) = ProcessType.ID_IMMEDIATE;
                            case 'DISABLED'
                                sn.procid(i,j) = ProcessType.ID_DISABLED;
                            case 'COX2'
                                sn.procid(i,j) = ProcessType.ID_COX2;
                        end
                    end
                end
            else
                sn.procid = [];
            end

            if ~isempty(java_sn.mu)
                sn.mu = cell(sn.nstations, 1);
                for i = 1:sn.nstations
                    sn.mu{i} = cell(1, sn.nclasses);
                    for j = 1:sn.nclasses
                        sn.mu{i}{j} = JLINE.jlinematrix_to_matrix(java_sn.mu.get(java_stations.get(i-1)).get(java_classes.get(j-1)));
                    end
                end
            else
                sn.mu = {};
            end

            if ~isempty(java_sn.phi)
                sn.phi = cell(sn.nstations, 1);
                for i = 1:sn.nstations
                    sn.phi{i} = cell(1, sn.nclasses);
                    for j = 1:sn.nclasses
                        sn.phi{i}{j} = JLINE.jlinematrix_to_matrix(java_sn.phi.get(java_stations.get(i-1)).get(java_classes.get(j-1)));
                    end
                end
            else
                sn.phi = {};
            end

            if ~isempty(java_sn.proc)
                sn.proc = cell(sn.nstations, 1);
                for i = 1:sn.nstations
                    sn.proc{i} = cell(1, sn.nclasses);
                    for j = 1:sn.nclasses
                        proc_i_j = java_sn.proc.get(java_stations.get(i-1)).get(java_classes.get(j-1));
                        sn.proc{i}{j} = cell(1, proc_i_j.size);
                        for k = 1:proc_i_j.size
                            sn.proc{i}{j}{k} = JLINE.jlinematrix_to_matrix(proc_i_j.get(uint32(k-1)));
                        end
                    end
                end
            else
                sn.proc = {};
            end

            if ~isempty(java_sn.pie)
                sn.pie = cell(sn.nstations, 1);
                for i = 1:sn.nstations
                    sn.pie{i} = cell(1, sn.nclasses);
                    for j = 1:sn.nclasses
                        sn.pie{i}{j} = JLINE.jlinematrix_to_matrix(java_sn.pie.get(java_stations.get(i-1)).get(java_classes.get(j-1)));
                    end
                end
            else
                sn.pie = {};
            end

            if ~isempty(java_sn.sched)
                sn.sched = cell(sn.nstations, 1);
                sn.schedid = zeros(sn.nstations, 1);
                for i = 1:sn.nstations
                    schedStrategy = java_sn.sched.get(java_stations.get(i-1));
                    switch schedStrategy.name.toCharArray'
                        case 'INF'
                            sn.sched{i} = SchedStrategy.INF;
                            sn.schedid(i) = SchedStrategy.ID_INF;
                        case 'FCFS'
                            sn.sched{i} = SchedStrategy.FCFS;
                            sn.schedid(i) = SchedStrategy.ID_FCFS;
                        case 'LCFS'
                            sn.sched{i} = SchedStrategy.LCFS;
                            sn.schedid(i) = SchedStrategy.ID_LCFS;
                        case 'LCFSPR'
                            sn.sched{i} = SchedStrategy.LCFSPR;
                            sn.schedid(i) = SchedStrategy.ID_LCFSPR;
                        case 'SIRO'
                            sn.sched{i} = SchedStrategy.SIRO;
                            sn.schedid(i) = SchedStrategy.ID_SIRO;
                        case 'SJF'
                            sn.sched{i} = SchedStrategy.SJF;
                            sn.schedid(i) = SchedStrategy.ID_SJF;
                        case 'LJF'
                            sn.sched{i} = SchedStrategy.LJF;
                            sn.schedid(i) = SchedStrategy.ID_LJF;
                        case 'PS'
                            sn.sched{i} = SchedStrategy.PS;
                            sn.schedid(i) = SchedStrategy.ID_PS;
                        case 'DPS'
                            sn.sched{i} = SchedStrategy.DPS;
                            sn.schedid(i) = SchedStrategy.ID_DPS;
                        case 'GPS'
                            sn.sched{i} = SchedStrategy.GPS;
                            sn.schedid(i) = SchedStrategy.ID_GPS;
                        case 'SEPT'
                            sn.sched{i} = SchedStrategy.SEPT;
                            sn.schedid(i) = SchedStrategy.ID_SEPT;
                        case 'LEPT'
                            sn.sched{i} = SchedStrategy.LEPT;
                            sn.schedid(i) = SchedStrategy.ID_LEPT;
                        case 'HOL'
                            sn.sched{i} = SchedStrategy.HOL;
                            sn.schedid(i) = SchedStrategy.ID_HOL;
                        case 'FORK'
                            sn.sched{i} = SchedStrategy.FORK;
                            sn.schedid(i) = SchedStrategy.ID_FORK;
                        case 'EXT'
                            sn.sched{i} = SchedStrategy.EXT;
                            sn.schedid(i) = SchedStrategy.ID_EXT;
                        case 'REF'
                            sn.sched{i} = SchedStrategy.REF;
                            sn.schedid(i) = SchedStrategy.ID_REF;
                    end
                end
            else
                sn.sched = {};
                sn.shcedid = [];
            end

            if ~isempty(java_sn.inchain)
                sn.inchain = cell(1, sn.nchains);
                for i = 1:sn.nchains
                    sn.inchain{1,i} = JLINE.jlinematrix_to_matrix(java_sn.inchain.get(uint32(i-1)))+1;
                end
            else
                sn.inchain = {};
            end

            if ~isempty(java_sn.visits)
                sn.visits = cell(sn.nchains, 1);
                for i = 1:sn.nchains
                    sn.visits{i,1} = JLINE.jlinematrix_to_matrix(java_sn.visits.get(uint32(i-1)));
                end
            else
                sn.visits = {};
            end

            if ~isempty(java_sn.nodevisits)
                sn.nodevisits = cell(1, sn.nchains);
                for i = 1:sn.nchains
                    sn.nodevisits{1,i} = JLINE.jlinematrix_to_matrix(java_sn.nodevisits.get(uint32(i-1)));
                end
            else
                sn.nodevisits = {};
            end

            if ~isempty(java_sn.droprule)
                sn.dropid = zeros(sn.nstations, sn.nclasses);
                for i = 1:sn.nstations
                    for j = 1:sn.nclasses
                        dropStrategy = java_sn.droprule.get(java_stations.get(i-1)).get(java_classes.get(j-1));
                        switch dropStrategy.name.toCharArray'
                            case 'WaitingQueue'
                                sn.dropid(i,j) = DropStrategy.ID_WAITQ;
                            case 'Drop'
                                sn.dropid(i,j) = DropStrategy.ID_DROP;
                            case 'BlockingAfterService'
                                sn.dropid(i,j) = DropStrategy.ID_BAS;
                        end
                    end
                end
            else
                sn.dropid = [];
            end

            if ~isempty(java_sn.nodeparam)
                sn.nodeparam = cell(sn.nnodes, 1);
                %Note that JLINE only support node parameters related to
                %Fork, Join, WWROBIN and RROBIN
                for i = 1:sn.nnodes
                    if java_sn.nodeparam.get(java_nodes.get(i-1)).isEmpty
                        sn.nodeparam{i} = [];
                    else
                        if ~isnan(java_sn.nodeparam.get(java_nodes.get(i-1)).fanOut)
                            sn.nodeparam{i}.fanOut = java_sn.nodeparam.get(java_nodes.get(i-1)).fanOut;
                        end
                        if ~isempty(java_sn.nodeparam.get(java_nodes.get(i-1)).joinStrategy)
                            if ~java_sn.nodeparam.get(java_nodes.get(i-1)).joinStrategy.isEmpty
                                sn.nodeparam{i}.joinStrategy = cell(1, sn.nclasses);
                                sn.nodeparam{i}.fanIn = cell(1, sn.nclasses);
                                for r = 1:sn.nclasses
                                    joinStrategy = java_sn.nodeparam.get(java_nodes.get(i-1)).joinStrategy.get(java_classes.get(r-1));
                                    switch joinStrategy.name.toCharArray'
                                        case 'STD'
                                            sn.nodeparam{i}.joinStrategy{r} = JoinStrategy.STD;
                                        case 'PARTIAL'
                                            sn.nodeparam{i}.joinStrategy{r} = JoinStrategy.PARTIAL;
                                    end
                                    sn.nodeparam{i}.fanIn{r} = java_sn.nodeparam.get(java_nodes.get(i-1)).fanIn.get(java_classes.get(r-1));
                                end
                            end
                        end

                        if ~isempty(java_sn.nodeparam.get(java_nodes.get(i-1)).weights)
                            for r = 1:sn.nclasses
                                sn.nodeparam{i}{r}.weights = JLINE.jlinematrix_to_matrix(java_sn.nodeparam.weights.get(java_classes.get(r-1)));
                            end
                        end

                        if ~isempty(java_sn.nodeparam.get(java_nodes.get(i-1)).outlinks)
                            for r = 1:sn.nclasses
                                sn.nodeparam{i}{r}.outlinks = JLINE.jlinematrix_to_matrix(java_sn.nodeparam.outlinks.get(java_classes.get(r-1)));
                            end
                        end
                    end
                end
            else
                sn.nodeparam = {};
            end

            if ~isempty(java_sn.sync)
                java_sync = java_sn.sync;
                sn.sync = cell(java_sync.size, 1);
                for i = 1:java_sync.size
                    java_sync_i = java_sync.get(uint32(i-1));
                    sn.sync{i,1} = struct('active',cell(1),'passive',cell(1));

                    java_active = java_sync_i.active.get(uint32(0));
                    java_passive = java_sync_i.passive.get(uint32(0));

                    %Currently assume that prob would always be a value
                    %instead of lambda function (No idea of how to convert
                    %Java lambda function to matlab lambda function)
                    switch java_active.getEvent.name.toCharArray'
                        case 'INIT'
                            sn.sync{i,1}.active{1} = Event(EventType.ID_INIT, java_active.getNodeIdx+1, java_active.getJobclassIdx+1, ...
                                java_active.getProb, JLINE.jlinematrix_to_matrix(java_active.getState), ...
                                java_active.getT, java_active.getJob);
                        case 'LOCAL'
                            sn.sync{i,1}.active{1} = Event(EventType.ID_LOCAL, java_active.getNodeIdx+1, java_active.getJobclassIdx+1, ...
                                java_active.getProb, JLINE.jlinematrix_to_matrix(java_active.getState), ...
                                java_active.getT, java_active.getJob);
                        case 'ARV'
                            sn.sync{i,1}.active{1} = Event(EventType.ID_ARV, java_active.getNodeIdx+1, java_active.getJobclassIdx+1, ...
                                java_active.getProb, JLINE.jlinematrix_to_matrix(java_active.getState), ...
                                java_active.getT, java_active.getJob);
                        case 'DEP'
                            sn.sync{i,1}.active{1} = Event(EventType.ID_DEP, java_active.getNodeIdx+1, java_active.getJobclassIdx+1, ...
                                java_active.getProb, JLINE.jlinematrix_to_matrix(java_active.getState), ...
                                java_active.getT, java_active.getJob);
                        case 'PHASE'
                            sn.sync{i,1}.active{1} = Event(EventType.ID_PHASE, java_active.getNodeIdx+1, java_active.getJobclassIdx+1, ...
                                java_active.getProb, JLINE.jlinematrix_to_matrix(java_active.getState), ...
                                java_active.getT, java_active.getJob);
                        case 'READ'
                            sn.sync{i,1}.active{1} = Event(EventType.ID_READ, java_active.getNodeIdx+1, java_active.getJobclassIdx+1, ...
                                java_active.getProb, JLINE.jlinematrix_to_matrix(java_active.getState), ...
                                java_active.getT, java_active.getJob);
                        case 'STAGE'
                            sn.sync{i,1}.active{1} = Event(EventType.ID_STAGE, java_active.getNodeIdx+1, java_active.getJobclassIdx+1, ...
                                java_active.getProb, JLINE.jlinematrix_to_matrix(java_active.getState), ...
                                java_active.getT, java_active.getJob);
                    end

                    switch java_passive.getEvent.name.toCharArray'
                        case 'INIT'
                            sn.sync{i,1}.passive{1} = Event(EventType.ID_INIT, java_passive.getNodeIdx+1, java_passive.getJobclassIdx+1, ...
                                java_passive.getProb, JLINE.jlinematrix_to_matrix(java_passive.getState), ...
                                java_passive.getT, java_passive.getJob);
                        case 'LOCAL'
                            sn.sync{i,1}.passive{1} = Event(EventType.ID_LOCAL, java_passive.getNodeIdx+1, java_passive.getJobclassIdx+1, ...
                                java_passive.getProb, JLINE.jlinematrix_to_matrix(java_passive.getState), ...
                                java_passive.getT, java_passive.getJob);
                        case 'ARV'
                            sn.sync{i,1}.passive{1} = Event(EventType.ID_ARV, java_passive.getNodeIdx+1, java_passive.getJobclassIdx+1, ...
                                java_passive.getProb, JLINE.jlinematrix_to_matrix(java_passive.getState), ...
                                java_passive.getT, java_passive.getJob);
                        case 'DEP'
                            sn.sync{i,1}.passive{1} = Event(EventType.ID_DEP, java_passive.getNodeIdx+1, java_passive.getJobclassIdx+1, ...
                                java_passive.getProb, JLINE.jlinematrix_to_matrix(java_passive.getState), ...
                                java_passive.getT, java_passive.getJob);
                        case 'PHASE'
                            sn.sync{i,1}.passive{1} = Event(EventType.ID_PHASE, java_passive.getNodeIdx+1, java_passive.getJobclassIdx+1, ...
                                java_passive.getProb, JLINE.jlinematrix_to_matrix(java_passive.getState), ...
                                java_passive.getT, java_passive.getJob);
                        case 'READ'
                            sn.sync{i,1}.passive{1} = Event(EventType.ID_READ, java_passive.getNodeIdx+1, java_passive.getJobclassIdx+1, ...
                                java_passive.getProb, JLINE.jlinematrix_to_matrix(java_passive.getState), ...
                                java_passive.getT, java_passive.getJob);
                        case 'STAGE'
                            sn.sync{i,1}.passive{1} = Event(EventType.ID_STAGE, java_passive.getNodeIdx+1, java_passive.getJobclassIdx+1, ...
                                java_passive.getProb, JLINE.jlinematrix_to_matrix(java_passive.getState), ...
                                java_passive.getT, java_passive.getJob);
                    end
                end
            else
                sn.sync = {};
            end
        end

        function [QN,UN,RN,WN,TN] = arrayListToResults(alist)
            n = alist.get(0).size;
            QN = zeros(n,1);
            for i=1:n
                QN(i) = alist.get(0).get(i-1);
            end
            UN = zeros(n,1);
            for i=1:n
                UN(i) = alist.get(1).get(i-1);
            end
            RN = zeros(n,1);
            for i=1:n
                RN(i) = alist.get(2).get(i-1);
            end
            WN = zeros(n,1);
            for i=1:n
                WN(i) = alist.get(3).get(i-1);
            end
            TN = zeros(n,1);
            for i=1:n
                TN(i) = alist.get(4).get(i-1);
            end
        end

        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            featSupported.setTrue({'Sink','Source',...
                'ClassSwitch','DelayStation','Queue',...
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

        function [ssa] = SolverSSA(network_object)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.SSA);
            ssa = jline.solvers.ssa.SolverSSA(network_object,solverOptions);
            ssa.disableResTime = true;
        end

        function [mam] = SolverMAM(network_object)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.MAM);
            mam = jline.solvers.mam.SolverMAM(network_object, solverOptions);
            %mam.compile(network_object);
        end

        function [nc] = SolverNC(network_object)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.NC);
            nc = jline.solvers.nc.SolverNC(network_object, solverOptions);
        end

        function [jmt] = SolverJMT(network_object)
            %solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.JMT);
            jmt = jline.solvers.jmt.SolverJMT(network_object);
            %jmt.compile(network_object);
        end

        function [ctmc] = SolverCTMC(network_object)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.CTMC);
            ctmc = jline.solvers.ctmc.SolverCTMC(network_object,solverOptions);
            ctmc.compile(network_object);
        end

        function [fluid] = SolverFluid(network_object, method)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.FLUID);
            solverOptions.verbose = solverOptions.verbose.SILENT;
            if nargin > 1
                solverOptions.method = method;
            end
            fluid = jline.solvers.fluid.SolverFluid(network_object, solverOptions);
        end

        function [mva] = SolverMVA(network_object, options)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.MVA);
            fn = fieldnames(options);
            fn2 = fieldnames(solverOptions);
            for f = 1:length(fn)
                found = 0;
                for j = 1:length(fn2)
                    if strcmp(fn{f}, fn2{j})
                        found = 1;
                        switch fn{f}
                            case 'config'
                                solverOptions.config.highvar = options.config.highvar;
                                solverOptions.config.multiserver = options.config.multiserver;
                                solverOptions.config.np_priority = options.config.np_priority;
                                solverOptions.config.fork_join = options.config.fork_join;
                            case 'verbose'
                                switch options.(fn{f})
                                    case 1
                                        solverOptions.verbose = solverOptions.verbose.SILENT;
                                    case 2
                                        solverOptions.verbose = solverOptions.verbose.STD;
                                    case 3
                                        solverOptions.verbose = solverOptions.verbose.DEBUG;
                                end
                            case 'init_sol'
                                solverOptions.(fn{f}) = JLINE.matrix_to_jlinematrix(options.init_sol);
                            case 'odesolvers'
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
            mva = jline.solvers.mva.SolverMVA(network_object, solverOptions);
        end

        function [ln] = SolverLN(layered_network_object)
            solverOptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.LN);
            ln = jline.solvers.ln.SolverLN(layered_network_object, solverOptions);
        end

    end
end
