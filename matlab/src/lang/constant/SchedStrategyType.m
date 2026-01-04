classdef (Sealed) SchedStrategyType
    % Enumeration of scheduling strategy types
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        NP = 0;        % Non-preemptive
        PR = 1;        % Preemptive resume
        PNR = 2;       % Preemptive non-resume
        NPPrio = 3;    % Non-preemptive priority
        PRPrio = 4;    % Preemptive resume priority
        PNRPrio = 5;   % Preemptive non-resume priority
    end

    methods (Static)
        function typeId = getTypeId(strategy)
            % TYPEID = GETTYPEID(STRATEGY)
            % Classifies the scheduling strategy type
            switch strategy
                case {SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, SchedStrategy.LPS}
                    typeId = SchedStrategyType.PR;
                case {SchedStrategy.FCFS, SchedStrategy.SIRO, ...
                      SchedStrategy.SEPT, SchedStrategy.LEPT, ...
                      SchedStrategy.SJF, SchedStrategy.INF, SchedStrategy.SETF}
                    typeId = SchedStrategyType.NP;
                case SchedStrategy.HOL
                    typeId = SchedStrategyType.NPPrio;
                case SchedStrategy.PSPRIO
                    typeId = SchedStrategyType.PRPrio;
                otherwise
                    line_error(mfilename, 'Unrecognized scheduling strategy type.');
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
                otherwise
                    line_error(mfilename, 'Unrecognized scheduling strategy type.');
            end
        end
    end

    methods (Access = private)
        function out = SchedStrategyType
            % Prevent instantiation
        end
    end
end
