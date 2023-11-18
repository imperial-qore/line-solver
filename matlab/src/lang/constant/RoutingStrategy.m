classdef (Sealed) RoutingStrategy
    % Enumeration of routing strategies
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    properties (Constant)
        ID_RAND = 0;
        ID_PROB = 1;
        ID_RROBIN = 2;
        ID_WRROBIN = 3;
        ID_JSQ = 4;
        ID_FIRING = 5;
        ID_KCHOICES = 6;
        ID_RL = 7;
        ID_DISABLED = -1;
        RAND = 'Random';
        RROBIN = 'RoundRobin';
        WRROBIN = 'WeightedRoundRobin';
        PROB = 'Probabilities';
        JSQ = 'JoinShortestQueue';
        DISABLED = 'Disabled';
        FIRING = 'Firing';
        KCHOICES = 'PowerKChoices';
        RL = 'ReinforcementLearning';
    end

    methods (Static, Access = public)
        function type = toType(text)
            % TYPE = TOTYPE(TEXT)

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
                case 'PowerKChoices'
                    type = RoutingStrategy.KCHOICES;
                case 'ReinforcementLearning'
                    type = RoutingStrategy.RL;
                case 'Disabled'
                    type = RoutingStrategy.DISABLED;
                otherwise
                    type = RoutingStrategy.ID_DISABLED;                    
            end
        end

        function id = toId(type)
            switch type
                case RoutingStrategy.RAND
                    id = RoutingStrategy.ID_RAND;
                case RoutingStrategy.PROB
                    id = RoutingStrategy.ID_PROB;
                case RoutingStrategy.RROBIN
                    id = RoutingStrategy.ID_RROBIN;
                case RoutingStrategy.WRROBIN
                    id = RoutingStrategy.ID_WRROBIN;
                case RoutingStrategy.FIRING
                    id = RoutingStrategy.ID_FIRING;
                case RoutingStrategy.JSQ
                    id = RoutingStrategy.ID_JSQ;
                case RoutingStrategy.KCHOICES
                    id = RoutingStrategy.ID_KCHOICES;
                case RoutingStrategy.RL
                    id = RoutingStrategy.ID_RL;
                case RoutingStrategy.DISABLED
                    id = RoutingStrategy.ID_DISABLED;
                otherwise
                    id = RoutingStrategy.ID_DISABLED;
            end
        end

        function feature = toFeature(type)
            % FEATURE = TOFEATURE(TYPE)
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
                    feature = 'RoutingStrategy_FIRING';
                case RoutingStrategy.JSQ
                    feature = 'RoutingStrategy_JSQ';
                case RoutingStrategy.KCHOICES
                    feature = 'RoutingStrategy_KCHOICES';
                case RoutingStrategy.RL
                    feature = 'RoutingStrategy_RL';
                case RoutingStrategy.DISABLED
                    feature = 'RoutingStrategy_DISABLED';
                otherwise % if unassigned, set it by default to Disabled
                    feature = 'RoutingStrategy_DISABLED';
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
                case RoutingStrategy.KCHOICES
                    text = 'KChoices';
                case RoutingStrategy.RL
                    text = 'RL';
                case RoutingStrategy.DISABLED
                    text = 'Disabled';
                otherwise % if unassigned, set it by default to Disabled
                    text = 'Disabled';
            end
        end
    end

end

