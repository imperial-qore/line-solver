classdef (Sealed) RoutingStrategy
    % RoutingStrategy Enumeration of job routing policies and load balancing strategies
    %
    % RoutingStrategy defines constants for routing policies that determine how
    % jobs are directed from one node to another in queueing networks. These
    % strategies control load distribution, traffic balancing, and path selection
    % throughout the network topology.
    %
    % @brief Comprehensive enumeration of job routing and load balancing strategies
    %
    % Key routing categories:
    % - Probabilistic: RAND, PROB (random and probability-based routing)
    % - Load balancing: JSQ, SQ (queue length-based decisions)
    % - Round-robin: RROBIN, WRROBIN (cyclic and weighted distribution)
    % - Advanced: RL, FIRING (learning and event-based routing)
    % - Control: DISABLED (no routing for specific classes)
    %
    % Common routing strategies:
    % - RAND: Random routing (equal probability to all destinations)
    % - PROB: Probabilistic routing (user-specified probabilities)
    % - RROBIN: Round-robin (cyclic distribution)
    % - JSQ: Join Shortest Queue (dynamic load balancing)
    % - SQ: Shortest queue of d, SQ(d) (formerly KCHOICES)
    % - DISABLED: No routing (class blocked at this node)
    %
    % RoutingStrategy is used in:
    % - Node output section configuration
    % - Network topology specification
    % - Load balancing implementation
    % - Traffic distribution control
    % - Router and dispatcher configuration
    %
    % Example:
    % @code
    % router.setRouting(jobClass, RoutingStrategy.PROB, [queue1, queue2], [0.7, 0.3]);
    % loadBalancer.setRouting(jobClass, RoutingStrategy.JSQ);
    % roundRobin.setRouting(jobClass, RoutingStrategy.RROBIN);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        RAND      = 0;
        PROB      = 1;
        RROBIN    = 2;
        WRROBIN   = 3;
        JSQ       = 4;
        FIRING    = 5;
        SQ        = 6;  % KCHOICES is now SQ: shortest queue of d, SQ(d)
        RL        = 7;
        DISABLED  = -1;
    end

    methods (Static, Access = public)

        function type = fromText(text)
            % TYPE = FROMTEXT(TEXT)
            switch text
                case 'Random'
                    type = RoutingStrategy.RAND;
                case 'Probabilities'
                    type = RoutingStrategy.PROB;
                case 'RoundRobin'
                    type = RoutingStrategy.RROBIN;
                case 'WeightedRoundRobin'
                    type = RoutingStrategy.WRROBIN;
                case 'JoinShortestQueue'
                    type = RoutingStrategy.JSQ;
                case 'Firing'
                    type = RoutingStrategy.FIRING;
                case 'ShortestQueueOfD'
                    type = RoutingStrategy.SQ;
                case 'ReinforcementLearning'
                    type = RoutingStrategy.RL;
                case 'Disabled'
                    type = RoutingStrategy.DISABLED;
                otherwise
                    line_error(mfilename, 'Unrecognized routing strategy string.');
            end
        end

        function id = toId(type)
            % ID = TOID(TYPE)
            switch type
                case RoutingStrategy.RAND
                    id = RoutingStrategy.RAND;
                case RoutingStrategy.PROB
                    id = RoutingStrategy.PROB;
                case RoutingStrategy.RROBIN
                    id = RoutingStrategy.RROBIN;
                case RoutingStrategy.WRROBIN
                    id = RoutingStrategy.WRROBIN;
                case RoutingStrategy.FIRING
                    id = RoutingStrategy.FIRING;
                case RoutingStrategy.JSQ
                    id = RoutingStrategy.JSQ;
                case RoutingStrategy.SQ
                    id = RoutingStrategy.SQ;
                case RoutingStrategy.RL
                    id = RoutingStrategy.RL;
                case RoutingStrategy.DISABLED
                    id = RoutingStrategy.DISABLED;
                otherwise
                    line_error(mfilename, 'Unrecognized routing strategy ID.');
            end
        end

        function feature = toFeature(type)
            % FEATURE = TOFEATURE(TYPE)
            if ischar(type) || isstring(type)
                % output.outputStrategy{r}{2} is heterogeneously typed: the
                % section constructors (Linkage, Forker, Firing) leave the
                % numeric id there, while Node.setRouting and Network.link
                % overwrite it with RoutingStrategy.toText(...). Accept both.
                % Without this, getUsedLangFeatures errored on every node whose
                % routing had been assigned; its enclosing TRY swallowed the
                % error, so no model ever reported a routing capability and the
                % solver support gate never saw one (state-dependent routing
                % passed silently to solvers that cannot represent it).
                type = RoutingStrategy.fromText(char(type));
            end
            switch type
                case RoutingStrategy.RAND
                    feature = 'RoutingStrategy_RAND';
                case RoutingStrategy.PROB
                    feature = 'RoutingStrategy_PROB';
                case RoutingStrategy.RROBIN
                    feature = 'RoutingStrategy_RROBIN';
                case RoutingStrategy.WRROBIN
                    feature = 'RoutingStrategy_WRROBIN';
                case RoutingStrategy.FIRING
                    % FIRING (Petri-net transition output) and DISABLED (class
                    % not routed at this node) are internal markers, not
                    % user-selectable routing capabilities, so they map to no
                    % registry name and are skipped by setUsedLangFeature. Same
                    % convention as jar RoutingStrategy.toFeature, which
                    % returns "" for both.
                    feature = '';
                case RoutingStrategy.JSQ
                    feature = 'RoutingStrategy_JSQ';
                case RoutingStrategy.SQ
                    feature = 'RoutingStrategy_SQ';
                case RoutingStrategy.RL
                    feature = 'RoutingStrategy_RL';
                case RoutingStrategy.DISABLED
                    feature = '';
                otherwise
                    line_error(mfilename, 'Unrecognized routing strategy feature.');
            end
        end

        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            switch type
                case RoutingStrategy.RAND
                    text = 'Random';
                case RoutingStrategy.PROB
                    text = 'Probabilities';
                case RoutingStrategy.RROBIN
                    text = 'RoundRobin';
                case RoutingStrategy.WRROBIN
                    text = 'WeightedRoundRobin';
                case RoutingStrategy.FIRING
                    text = 'Firing';
                case RoutingStrategy.JSQ
                    text = 'JoinShortestQueue';
                case RoutingStrategy.SQ
                    text = 'ShortestQueueOfD';
                case RoutingStrategy.RL
                    text = 'ReinforcementLearning';
                case RoutingStrategy.DISABLED
                    text = 'Disabled';
                otherwise
                    line_error(mfilename, 'Unrecognized routing strategy type.');
            end
        end
    end
end
