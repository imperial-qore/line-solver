classdef (Sealed) ProcessType
    % Enumeration of process ts
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        EXP = 'exp'; % renewal
        ERLANG = 'erlang';% renewal
        HYPEREXP = 'hyperexp';% renewal
        PH = 'ph'; % renewal
        APH = 'aph';% renewal
        MAP = 'map';% renewal
        UNIFORM = 'uniform';% renewal
        DET = 'det';% renewal
        COXIAN = 'coxian';% renewal
        GAMMA = 'gamma';% renewal
        PARETO = 'pareto';% renewal
        WEIBULL = 'weibull';% renewal
        LOGNORMAL = 'lognormal';% renewal
        MMPP2 = 'mmpp2';
        REPLAYER = 'trace';
        TRACE = 'trace';
        IMMEDIATE = 'immediate';
        DISABLED = 'disabled';
        COX2 = 'cox2';
        
        ID_EXP = 0;
        ID_ERLANG = 1;
        ID_HYPEREXP = 2;
        ID_PH = 3;
        ID_APH = 4;
        ID_MAP = 5;
        ID_UNIFORM = 6;
        ID_DET = 7;
        ID_COXIAN = 8;
        ID_GAMMA = 9;
        ID_PARETO = 10;
        ID_MMPP2 = 11;
        ID_REPLAYER = 12;
        ID_TRACE = 12;
        ID_IMMEDIATE = 13;
        ID_DISABLED = 14;
        ID_COX2 = 15;
        ID_WEIBULL = 16;
        ID_LOGNORMAL = 17;
    end
    
    methods (Static)        
        function t = fromId(id)
            % ID = TOID(TYPE)
            switch id
                case ProcessType.ID_EXP
                    t = ProcessType.EXP;
                case ProcessType.ID_ERLANG
                    t = ProcessType.ERLANG;
                case ProcessType.ID_HYPEREXP
                    t = ProcessType.HYPEREXP;
                case ProcessType.ID_PH
                    t = ProcessType.PH;
                case ProcessType.ID_APH
                    t = ProcessType.APH;
                case ProcessType.ID_MAP
                    t = ProcessType.MAP;
                case ProcessType.ID_UNIFORM
                    t = ProcessType.UNIFORM;
                case ProcessType.ID_DET
                    t = ProcessType.DET;
                case ProcessType.ID_COXIAN
                    t = ProcessType.COXIAN;
                case ProcessType.ID_GAMMA
                    t = ProcessType.GAMMA;
                case ProcessType.ID_PARETO
                    t = ProcessType.PARETO;
                case ProcessType.ID_MMPP2
                    t = ProcessType.MMPP2;
                case ProcessType.ID_TRACE
                    t = ProcessType.TRACE;
                case ProcessType.ID_IMMEDIATE
                    t = ProcessType.IMMEDIATE;
                case ProcessType.ID_DISABLED
                    t = ProcessType.DISABLED;
                case ProcessType.ID_WEIBULL
                    t = ProcessType.WEIBULL;
                case ProcessType.ID_LOGNORMAL
                    t = ProcessType.LOGNORMAL;
            end
        end
        
        function id = toId(t)
            % ID = TOID(TYPE)            
            switch t
                case ProcessType.EXP
                    id = 0;
                case ProcessType.ERLANG
                    id = 1;
                case ProcessType.HYPEREXP
                    id = 2;
                case ProcessType.PH
                    id = 3;
                case ProcessType.APH
                    id = 4;
                case ProcessType.MAP
                    id = 5;
                case ProcessType.UNIFORM
                    id = 6;
                case ProcessType.DET
                    id = 7;
                case ProcessType.COXIAN
                    id = 8;
                case ProcessType.GAMMA
                    id = 9;
                case ProcessType.PARETO
                    id = 10;
                case ProcessType.MMPP2
                    id = 11;
                case {ProcessType.REPLAYER, ProcessType.TRACE}
                    id = 12;
                case ProcessType.IMMEDIATE
                    id = 13;
                case ProcessType.DISABLED
                    id = 14;
                case ProcessType.COX2
                    id = 15;
                case ProcessType.WEIBULL
                    id = 16;
                case ProcessType.LOGNORMAL
                    id = 17;
            end
        end
        
        function t = fromText(text)
            % TIMMEDIATE = TOID(TYPE)
            switch text
                case 'Exp'
                    t = ProcessType.EXP;
                case 'Erlang'
                    t = ProcessType.ERLANG;
                case 'HyperExp'
                    t = ProcessType.HYPEREXP;
                case 'PH'
                    t = ProcessType.PH;
                case 'APH'
                    t = ProcessType.APH;
                case 'MAP'
                    t = ProcessType.MAP;
                case 'Uniform'
                    t = ProcessType.UNIFORM;
                case 'Det'
                    t = ProcessType.DET;
                case 'Coxian'
                    t = ProcessType.COXIAN;
                case 'Gamma'
                    t = ProcessType.GAMMA;
                case 'Pareto'
                    t = ProcessType.PARETO;
                case 'MMPP2'
                    t = ProcessType.MMPP2;
                case 'Replayer'
                    t = ProcessType.TRACE;
                case 'Immediate'
                    t = ProcessType.IMMEDIATE;
                case 'Disabled'
                    t = ProcessType.DISABLED;
                case 'Cox2'
                    t = ProcessType.COX2;
                case 'Weibull'
                    t = ProcessType.WEIBULL;
                case 'Lognormal'
                    t = ProcessType.LOGNORMAL;
            end
        end
        
        function text = toText(t)
            % TEXT = TOTEXT(TYPE)
            switch t
                case ProcessType.EXP
                    text = 'Exp';
                case ProcessType.ERLANG
                    text = 'Erlang';
                case ProcessType.HYPEREXP
                    text = 'HyperExp';
                case ProcessType.PH
                    text = 'PH';
                case ProcessType.APH
                    text = 'APH';
                case ProcessType.MAP
                    text = 'MAP';
                case ProcessType.UNIFORM
                    text = 'Uniform';
                case ProcessType.DET
                    text = 'Det';
                case ProcessType.COXIAN
                    text = 'Coxian';
                case ProcessType.GAMMA
                    text = 'Gamma';
                case ProcessType.PARETO
                    text = 'Pareto';
                case ProcessType.MMPP2
                    text = 'MMPP2';
                case {ProcessType.REPLAYER, ProcessType.TRACE}
                    text = 'Replayer';
                case ProcessType.IMMEDIATE
                    text = 'Immediate';
                case ProcessType.DISABLED
                    text = 'Disabled';
                case ProcessType.COX2
                    text = 'Cox2';
                case ProcessType.WEIBULL
                    text = 'Weibull';
                case ProcessType.LOGNORMAL
                    text = 'Lognormal';
            end
            
        end
        
    end
end
