classdef Cluster < handle
    % CLUSTER Builder for cluster models with comparison helpers.
    %
    % Mirrors the Java jline.gen.Cluster and the Python
    % line_solver.gen.Cluster builders.
    %
    % Example:
    %   cluster = Cluster('numStations', 4, 'arrivalRate', 1.0, 'serviceRate', 0.4);
    %   cluster.setDispatching(RoutingStrategy.RAND).setScheduling(SchedStrategy.PS);
    %   model = cluster.build();
    %   table = SolverMVA(model).getAvgTable();
    %
    % The cluster is open by default; setClosed makes every class closed and
    % setMixed combines open and closed classes, open ones first.

    properties
        numStations (1, 1) double = 2
        arrivalRates double = 1.0
        serviceRates double = 1.0   % (M, R) matrix of *rates* (1/mean service time)
        stationCounts double          % per-server multiplicity
        scheduling = SchedStrategy.PS
        dispatching = RoutingStrategy.RAND
        closed (1, 1) logical = false
        mixed (1, 1) logical = false  % open classes first, then closed ones
        population double = []        % per-class population (closed only)
        thinkTimes double = []        % per-class think time (closed only)
        dispatchProbs = []            % PROB: rows = classes (1 broadcasts), cols = servers
        dispatchWeights = []          % WRROBIN: per-server integer weights
        sqD = []                % SQ: d, number of sampled destinations
        arrivalScvs = []              % per-class arrival SCVs (1 = exponential, default)
        serviceScvs = []              % (M, R) service SCVs (1 = exponential, default)
    end

    methods
        function obj = Cluster()
            % CLUSTER Default cluster: 2 servers, single open class with arrival
            % rate 1.0, service rate 1.0 at each server, PS scheduling, RAND
            % dispatching. Configure further via the chainable set* methods.
            obj.numStations   = 2;
            obj.arrivalRates = 1.0;
            obj.serviceRates = ones(2, 1);
            obj.stationCounts = ones(2, 1);
            obj.scheduling   = SchedStrategy.PS;
            obj.dispatching  = RoutingStrategy.RAND;
            obj.closed       = false;
            obj.mixed        = false;
        end

        function [R, Ro, Rc] = numClasses(obj)
            % [R, RO, RC] = NUMCLASSES()  Total, open and closed class counts.
            if obj.mixed
                Ro = numel(obj.arrivalRates);
                Rc = numel(obj.population);
            elseif obj.closed
                Ro = 0;
                Rc = numel(obj.population);
            else
                Ro = numel(obj.arrivalRates);
                Rc = 0;
            end
            R = Ro + Rc;
        end

        function obj = setNumStations(obj, M)
            % SETNUMSTATIONS  Number of parallel server queues. Replicates the
            % current single-class service rate across all servers and resets
            % the per-server multiplicity to all-1.
            if M <= 0
                error('numStations must be positive');
            end
            if isempty(obj.serviceRates)
                sample = 1.0;
            else
                sample = obj.serviceRates(1, 1);
            end
            obj.numStations   = M;
            obj.serviceRates = repmat(sample, M, 1);
            obj.stationCounts = ones(M, 1);
        end

        function obj = setArrivalRate(obj, lambda)
            % SETARRIVALRATE  Single-class arrival rate.
            if ~isscalar(lambda) || lambda <= 0
                error('arrival rate must be a positive scalar');
            end
            obj.arrivalRates = lambda;
        end

        function obj = setArrivalRates(obj, lambdas)
            % SETARRIVALRATES  Per-class arrival rates (length R).
            if any(lambdas <= 0)
                error('arrival rates must be positive');
            end
            obj.arrivalRates = lambdas(:)';
        end

        function obj = setServiceRate(obj, mu)
            % SETSERVICERATE  Single service rate, broadcast across all servers.
            if ~isscalar(mu) || mu <= 0
                error('service rate must be a positive scalar');
            end
            obj.serviceRates = repmat(mu, obj.numStations, 1);
        end

        function obj = setServiceRates(obj, rates)
            % SETSERVICERATES  Per-(server, class) service rates as a (M x R) matrix.
            if size(rates, 1) ~= obj.numStations
                error('serviceRates outer dim must equal numStations');
            end
            if any(rates(:) <= 0)
                error('service rates must be positive');
            end
            obj.serviceRates = rates;
        end

        function obj = setDispatching(obj, dispatching)
            obj.dispatching = dispatching;
        end

        function obj = setScheduling(obj, scheduling)
            obj.scheduling = scheduling;
        end

        function obj = setStationServers(obj, counts)
            if numel(counts) ~= obj.numStations
                error('counts length must equal numStations');
            end
            obj.stationCounts = counts(:);
        end

        function obj = setProbabilities(obj, probs)
            % SETPROBABILITIES  PROB dispatching with per-server probabilities.
            %
            % If PROBS is a row vector of length numStations it is broadcast to
            % every class; otherwise it must be a (R x numStations) matrix.
            if isvector(probs)
                if numel(probs) ~= obj.numStations
                    error('probs length must equal numStations');
                end
                obj.dispatchProbs = probs(:)';   % single row, broadcast at build
            else
                if size(probs, 2) ~= obj.numStations
                    error('probs must have numStations columns');
                end
                obj.dispatchProbs = probs;
            end
            obj.dispatching = RoutingStrategy.PROB;
        end

        function obj = setWeights(obj, weights)
            % SETWEIGHTS  WRROBIN dispatching with per-server integer weights.
            if numel(weights) ~= obj.numStations
                error('weights length must equal numStations');
            end
            if any(abs(weights - round(weights)) > 1e-9)
                error('WRROBIN weights must be integers');
            end
            obj.dispatchWeights = weights(:)';
            obj.dispatching = RoutingStrategy.WRROBIN;
        end

        function obj = setArrivalSCV(obj, scv)
            % SETARRIVALSCV  SCV of the arrival process (open clusters only).
            %
            % Either a scalar (broadcast to all classes) or a row vector of
            % length R. SCV != 1 swaps the per-class Exp distribution for
            % APH.fitMeanAndSCV(1/rate, scv).
            R = numel(obj.arrivalRates);
            if isscalar(scv)
                obj.arrivalScvs = scv * ones(1, R);
            else
                if numel(scv) ~= R
                    error('arrivalSCV length must equal number of classes');
                end
                obj.arrivalScvs = scv(:)';
            end
            if any(obj.arrivalScvs <= 0)
                error('SCV must be positive');
            end
        end

        function obj = setServiceSCV(obj, scv)
            % SETSERVICESCV  SCV of the per-server service distribution.
            %
            % Either a scalar (broadcast), a length-M vector (per-server,
            % broadcast across classes), or a (M x R) matrix.
            R = obj.numClasses();
            M = obj.numStations;
            if isscalar(scv)
                obj.serviceScvs = scv * ones(M, R);
            elseif isvector(scv) && numel(scv) == M
                obj.serviceScvs = repmat(scv(:), 1, R);
            elseif size(scv, 1) == M && size(scv, 2) == R
                obj.serviceScvs = scv;
            else
                error('serviceSCV must be a scalar, length-M vector, or (M x R) matrix');
            end
            if any(obj.serviceScvs(:) <= 0)
                error('SCV must be positive');
            end
        end

        function obj = setSQ(obj, d)
            % SETSQ  SQ(d) dispatching: shortest of d sampled servers.
            if d < 1 || d ~= round(d)
                error('d must be a positive integer');
            end
            obj.sqD = d;
            obj.dispatching = RoutingStrategy.SQ;
        end

        function obj = setClosed(obj, population, thinkTime)
            obj.closed = true;
            obj.mixed = false;
            if isscalar(population)
                obj.population = population;
                obj.thinkTimes = thinkTime;
            else
                if numel(population) ~= numel(thinkTime)
                    error('population and thinkTime must have the same length');
                end
                obj.population = population(:)';
                obj.thinkTimes = thinkTime(:)';
            end
        end

        function obj = setMixed(obj, arrivalRates, population, thinkTime)
            % SETMIXED  Mixed cluster: open classes coexist with closed ones.
            %
            % ARRIVALRATES holds the per-class arrival rates of the open
            % classes, POPULATION and THINKTIME the per-class population and
            % think time of the closed classes. Classes are ordered open first,
            % so serviceRates/serviceSCV matrices have numel(arrivalRates) +
            % numel(population) columns.
            if isempty(arrivalRates) || isempty(population)
                error('a mixed cluster needs at least one open and one closed class');
            end
            if any(arrivalRates <= 0)
                error('arrival rates must be positive');
            end
            if numel(population) ~= numel(thinkTime)
                error('population and thinkTime must have the same length');
            end
            obj.arrivalRates = arrivalRates(:)';
            obj.population = population(:)';
            obj.thinkTimes = thinkTime(:)';
            obj.closed = false;
            obj.mixed = true;
        end

        function model = build(obj)
            % BUILD  Construct the configured Network model.
            M = obj.numStations;
            R = obj.numClasses();

            strategy = cell(M, 1);
            for i = 1:M
                strategy{i} = obj.scheduling;
            end

            % Replicate single-rate vector to all classes if needed.
            sr = obj.serviceRates;
            if size(sr, 2) == 1 && R > 1
                sr = repmat(sr, 1, R);
            end
            if any(sr(:) <= 0)
                error('service rates must be positive');
            end
            D = 1.0 ./ sr;

            % Strategies that need per-destination parameters cannot be passed
            % to the static factory (which calls setRouting(class, strategy)
            % with no extras). Build with RAND and reapply in post-processing.
            needsPostProcess = ...
                (obj.dispatching == RoutingStrategy.PROB && ~isempty(obj.dispatchProbs)) ...
                || (obj.dispatching == RoutingStrategy.WRROBIN && ~isempty(obj.dispatchWeights)) ...
                || (obj.dispatching == RoutingStrategy.SQ && ~isempty(obj.sqD));
            if needsPostProcess
                factoryDispatch = RoutingStrategy.RAND;
            else
                factoryDispatch = obj.dispatching;
            end

            if obj.mixed
                model = MNetwork.clusterMixed(obj.arrivalRates, obj.population, ...
                    obj.thinkTimes, D, strategy, obj.stationCounts, factoryDispatch);
            elseif obj.closed
                model = MNetwork.clusterClosed(obj.population, obj.thinkTimes, ...
                    D, strategy, obj.stationCounts, factoryDispatch);
            else
                model = MNetwork.cluster(obj.arrivalRates, D, strategy, ...
                    obj.stationCounts, factoryDispatch);
            end

            obj.applyDispatcherConfig(model, R);
            obj.applyDistributionScvs(model, R);
        end

        function applyDistributionScvs(obj, model, R)
            jobclasses = model.classes;
            % Arrival SCVs (open classes only, which come first in the order).
            if ~obj.closed && ~isempty(obj.arrivalScvs)
                src = model.getNodeByName('Source');
                for r = 1:numel(obj.arrivalRates)
                    scv = obj.arrivalScvs(r);
                    if scv ~= 1
                        src.setArrival(jobclasses{r}, ...
                            APH.fitMeanAndSCV(1.0 / obj.arrivalRates(r), scv));
                    end
                end
            end
            % Service SCVs.
            if ~isempty(obj.serviceScvs)
                sr = obj.serviceRates;
                if size(sr, 2) == 1 && R > 1
                    sr = repmat(sr, 1, R);
                end
                for i = 1:obj.numStations
                    server = model.getNodeByName(['Station', num2str(i)]);
                    for r = 1:R
                        scv = obj.serviceScvs(i, r);
                        if scv ~= 1
                            server.setService(jobclasses{r}, ...
                                APH.fitMeanAndSCV(1.0 / sr(i, r), scv));
                        end
                    end
                end
            end
        end

        function applyDispatcherConfig(obj, model, R)
            if obj.dispatching == RoutingStrategy.PROB && ~isempty(obj.dispatchProbs)
                dispatcher = model.getNodeByName('Dispatcher');
                jobclasses = model.classes;
                for r = 1:R
                    if size(obj.dispatchProbs, 1) == 1
                        probsRow = obj.dispatchProbs;
                    else
                        probsRow = obj.dispatchProbs(r, :);
                    end
                    for i = 1:obj.numStations
                        server = model.getNodeByName(['Station', num2str(i)]);
                        dispatcher.setProbRouting(jobclasses{r}, server, probsRow(i));
                    end
                end
            elseif obj.dispatching == RoutingStrategy.WRROBIN && ~isempty(obj.dispatchWeights)
                dispatcher = model.getNodeByName('Dispatcher');
                jobclasses = model.classes;
                for r = 1:R
                    for i = 1:obj.numStations
                        server = model.getNodeByName(['Station', num2str(i)]);
                        dispatcher.setRouting(jobclasses{r}, RoutingStrategy.WRROBIN, ...
                            server, obj.dispatchWeights(i));
                    end
                end
            elseif obj.dispatching == RoutingStrategy.SQ && ~isempty(obj.sqD)
                dispatcher = model.getNodeByName('Dispatcher');
                jobclasses = model.classes;
                for r = 1:R
                    dispatcher.setRouting(jobclasses{r}, RoutingStrategy.SQ, obj.sqD);
                end
            end
        end

        function out = compareDispatching(obj, solverFcn, policies)
            % OUT = COMPAREDISPATCHING(SOLVERFCN, POLICIES) returns a dictionary
            % from each policy in POLICIES to the AvgTable produced by SOLVERFCN(model).
            out = configureDictionary('string', 'cell');
            saved = obj.dispatching;
            cleaner = onCleanup(@() obj.restoreDispatching(saved));
            for k = 1:numel(policies)
                obj.dispatching = policies(k);
                out{string(policies(k))} = solverFcn(obj.build());
            end
        end

        function out = compareScheduling(obj, solverFcn, disciplines)
            out = configureDictionary('string', 'cell');
            saved = obj.scheduling;
            cleaner = onCleanup(@() obj.restoreScheduling(saved));
            for k = 1:numel(disciplines)
                obj.scheduling = disciplines(k);
                out{string(disciplines(k))} = solverFcn(obj.build());
            end
        end

        function out = sweepArrivalRate(obj, rates, solverFcn)
            if obj.closed
                error('sweepArrivalRate is only defined for clusters with open classes');
            end
            if numel(obj.arrivalRates) ~= 1
                error('sweepArrivalRate requires a single open class');
            end
            out = configureDictionary('double', 'cell');
            saved = obj.arrivalRates;
            cleaner = onCleanup(@() obj.restoreArrivalRates(saved));
            for r = rates(:)'
                obj.arrivalRates = r;
                out{r} = solverFcn(obj.build());
            end
        end

        function out = sweepNumStations(obj, counts, solverFcn)
            out = configureDictionary('double', 'cell');
            savedM = obj.numStations;
            savedRates = obj.serviceRates;
            savedCounts = obj.stationCounts;
            cleaner = onCleanup(@() obj.restoreServers(savedM, savedRates, savedCounts));
            perClass = savedRates(1, :);
            for m = counts(:)'
                obj.numStations = m;
                obj.serviceRates = repmat(perClass, m, 1);
                obj.stationCounts = ones(m, 1);
                out{m} = solverFcn(obj.build());
            end
        end
    end

    methods (Access = private)
        function restoreDispatching(obj, saved)
            obj.dispatching = saved;
        end
        function restoreScheduling(obj, saved)
            obj.scheduling = saved;
        end
        function restoreArrivalRates(obj, saved)
            obj.arrivalRates = saved;
        end
        function restoreServers(obj, M, rates, counts)
            obj.numStations = M;
            obj.serviceRates = rates;
            obj.stationCounts = counts;
        end
    end
end
