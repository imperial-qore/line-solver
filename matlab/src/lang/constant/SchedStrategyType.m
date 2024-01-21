classdef (Sealed) SchedStrategyType
    % Enumeration of scheduling strategy types
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        ID_NP = 0;
        ID_PR = 1;
        ID_PNR = 2;
        ID_NPPrio = 3;
        ID_PRPrio = 4;
        ID_PNRPrio = 5;
        
        NP = 'NPR'; % non-preemptive
        PR = 'PR'; % preemptive resume
        PNR = 'PNR'; % preemptive non-resume
        NPPrio = 'NPPrio'; % non-preemptive priority
        PRPrio = 'PRPrio'; % preemptive resume priority
        PNRPrio = 'PNRPrio'; % preemptive non-resume priority
    end
    
    methods (Access = private)
        %private so that it cannot be instatiated.
        function out = SchedStrategyType
            % OUT = SCHEDSTRATEGYTYPE
            
        end
    end
    
    methods (Static)
        function typeId = getTypeId(strategy)
            % TYPEID = GETTYPEID(STRATEGY)
            % Classifies the scheduling strategy type
            switch strategy
                case {SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS}
                    typeId = SchedStrategyType.ID_PR;
                case {SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.ID_SEPT, SchedStrategy.ID_LEPT, SchedStrategy.ID_SJF, SchedStrategy.ID_INF}
                    typeId = SchedStrategyType.ID_NP;
                case SchedStrategy.HOL
                    typeId = SchedStrategyType.ID_NPPrio;
                otherwise
                    line_error(mfilename,'Unrecognized scheduling strategy type.');
            end
        end
        
        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            switch type
                case SchedStrategyType.NP
                    text = 'NonPreemptive';
                case SchedStrategyType.PR
                    text = 'PreemptiveResume';
                case SchedStrategyType.PNR
                    text = 'PreemptiveNonResume';
                case SchedStrategyType.NPPrio
                    text = 'NonPreemptivePriority';
                case SchedStrategyType.PRPrio
                    text = 'PreemptiveResumePriority';
                case SchedStrategyType.PNRPrio
                    text = 'PreemptiveNonResumePriority';
            end
        end
    end
end
