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
            % Classifies the scheduling strategy type.
            %
            % This IS the classification stored in schedPolicy: Queue.m calls
            % it and keeps no table of its own, and the JAR
            % SchedStrategyType.getTypeId and the native python get_type_id are
            % the same table. The three were written out separately until
            % 2026-09-17 and disagreed on eleven disciplines. Nothing in any
            % codebase READS schedPolicy, so the disagreement was never a wrong
            % number, only a wrong label -- which is exactly why it survived.
            %
            % The preemption axis follows each discipline's own definition in
            % SchedStrategy.m. PI is "preemptive independent", i.e. the
            % preempted job RESTARTS, so it is non-resume and not resume;
            % SRPT, PSJF, FB, LRPT, FSP and EDF are preemptive there, against
            % SJF, LJF, SEPT, LEPT, EDD and SETF, which rank the same jobs
            % without preempting. A discipline that also carries a priority
            % order is reported by its Prio type rather than by the plain one.
            switch SchedStrategy.toId(strategy)
                case {SchedStrategy.INF, SchedStrategy.FCFS, SchedStrategy.LCFS, ...
                      SchedStrategy.SIRO, SchedStrategy.SJF, SchedStrategy.LJF, ...
                      SchedStrategy.SEPT, SchedStrategy.LEPT, SchedStrategy.EDD, ...
                      SchedStrategy.SETF, SchedStrategy.SET, ...
                      SchedStrategy.PAS, SchedStrategy.OI, SchedStrategy.POLLING}
                    typeId = SchedStrategyType.NP;
                case {SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, ...
                      SchedStrategy.LPS, SchedStrategy.LCFSPR, SchedStrategy.FCFSPR, ...
                      SchedStrategy.EDF, SchedStrategy.SRPT, SchedStrategy.PSJF, ...
                      SchedStrategy.FB, SchedStrategy.LAS, SchedStrategy.LRPT, ...
                      SchedStrategy.FSP}
                    typeId = SchedStrategyType.PR;
                case {SchedStrategy.FCFSPI, SchedStrategy.LCFSPI}
                    typeId = SchedStrategyType.PNR;
                case {SchedStrategy.HOL, SchedStrategy.FCFSPRIO, SchedStrategy.LCFSPRIO}
                    typeId = SchedStrategyType.NPPrio;
                case {SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, ...
                      SchedStrategy.LCFSPRPRIO, SchedStrategy.FCFSPRPRIO, ...
                      SchedStrategy.SRPTPRIO}
                    typeId = SchedStrategyType.PRPrio;
                case {SchedStrategy.FCFSPIPRIO, SchedStrategy.LCFSPIPRIO}
                    typeId = SchedStrategyType.PNRPrio;
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
