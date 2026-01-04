classdef (Sealed) SchedStrategy
    % SchedStrategy Enumeration of queueing scheduling disciplines and strategies
    %
    % SchedStrategy defines constants for all supported scheduling disciplines
    % in LINE queueing systems. These strategies determine the order in which
    % jobs are selected for service from queues, affecting system performance
    % and fairness characteristics.
    %
    % @brief Comprehensive enumeration of queueing scheduling disciplines
    %
    % Key scheduling categories:
    % - Order-based: FCFS, LCFS, SIRO (service order policies)
    % - Size-based: SJF, LJF, SEPT, LEPT (job length policies)
    % - Sharing: PS, DPS, GPS (processor sharing variants)
    % - Priority: HOL, PSPRIO, DPSPRIO, GPSPRIO (priority disciplines)
    % - Special: INF, FORK, POLLING, EXT, REF (specialized strategies)
    %
    % Common scheduling strategies:
    % - FCFS: First-Come-First-Served (FIFO)
    % - LCFS: Last-Come-First-Served (LIFO/Stack)
    % - PS: Processor Sharing (round-robin with infinitesimal time slices)
    % - SJF: Shortest Job First (preemptive shortest remaining time)
    % - HOL: Head-of-Line priority (non-preemptive priority)
    % - INF: Infinite server (delay station, no queueing)
    %
    % SchedStrategy is used in:
    % - Queue node configuration
    % - Service discipline specification
    % - Performance analysis parameterization
    % - Model validation and compatibility checking
    % - Solver method selection
    %
    % Example:
    % @code
    % queue1 = Queue(model, 'Server1', SchedStrategy.FCFS);
    % queue2 = Queue(model, 'Server2', SchedStrategy.PS);
    % delay = Queue(model, 'ThinkTime', SchedStrategy.INF);
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        INF      = 0;
        FCFS     = 1;
        LCFS     = 2;
        SIRO     = 3;
        SJF      = 4;
        LJF      = 5;
        PS       = 6;
        DPS      = 7;
        GPS      = 8;
        SEPT     = 9;
        LEPT     = 10;
        HOL      = 11;
        FORK     = 12;
        EXT      = 13;
        REF      = 14;
        LCFSPR   = 15;
        POLLING  = 16;
        PSPRIO   = 17;
        DPSPRIO  = 18;
        GPSPRIO  = 19;
        LCFSPI   = 20;
        LCFSPRIO = 21;
        LCFSPRPRIO = 22;
        LCFSPIPRIO = 23;
        FCFSPRIO = 11;  % Alias for HOL
        FCFSPR   = 24;
        FCFSPI   = 25;
        FCFSPRPRIO = 26;
        FCFSPIPRIO = 27;
        SRPT     = 28;
        SRPTPRIO = 29;
        EDD      = 30;  % Earliest Due Date (non-preemptive)
        EDF      = 31;  % Earliest Deadline First (preemptive)
        LPS      = 32;  % Least Progress Scheduling
        PSJF     = 33;  % Preemptive Shortest Job First (priority by original size)
        FB       = 34;  % Feedback / Least Attained Service (priority by age)
        LAS      = 34;  % Alias for FB (Least Attained Service)
        LRPT     = 35;  % Longest Remaining Processing Time
        SETF     = 36;  % Shortest Elapsed Time First (non-preemptive FB)
        SET      = 36;  % Alias for SETF
    end

    methods (Static)
        function id = toId(type)
            % ID = TOID(TYPE)
            if isnumeric(type)
                id = type;
            else
                line_error(mfilename, 'Invalid scheduling strategy type.');
            end
        end

        function type = fromId(id)
            % TYPE = FROMID(ID)
            if ismember(id, [SchedStrategy.INF, SchedStrategy.FCFS, SchedStrategy.LCFS, ...
                             SchedStrategy.SIRO, SchedStrategy.SJF, SchedStrategy.LJF, ...
                             SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, ...
                             SchedStrategy.SEPT, SchedStrategy.LEPT, SchedStrategy.SRPT, SchedStrategy.SRPTPRIO, SchedStrategy.HOL, ...
                             SchedStrategy.FORK, SchedStrategy.EXT, SchedStrategy.REF, ...
                             SchedStrategy.LCFSPR, SchedStrategy.POLLING, ...
                             SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, ...
                             SchedStrategy.LCFSPI, SchedStrategy.LCFSPRIO, SchedStrategy.LCFSPRPRIO, SchedStrategy.LCFSPIPRIO, ...
                             SchedStrategy.FCFSPRIO, SchedStrategy.FCFSPR, SchedStrategy.FCFSPI, SchedStrategy.FCFSPRPRIO, SchedStrategy.FCFSPIPRIO, ...
                             SchedStrategy.EDD, SchedStrategy.EDF, SchedStrategy.LPS, ...
                            SchedStrategy.PSJF, SchedStrategy.FB, SchedStrategy.LRPT, SchedStrategy.SETF])
                type = id;
            else
                line_error(mfilename, 'Unrecognized scheduling strategy ID.');
            end
        end

        function text = toText(type)
            % TEXT = TOTEXT(TYPE)
            switch type
                case SchedStrategy.INF
                    text = 'inf';
                case SchedStrategy.FCFS
                    text = 'fcfs';
                case SchedStrategy.FCFSPR
                    text = 'fcfspr';
                case SchedStrategy.FCFSPI
                    text = 'fcfspi';
                case SchedStrategy.LCFS
                    text = 'lcfs';
                case SchedStrategy.SIRO
                    text = 'siro';
                case SchedStrategy.SJF
                    text = 'sjf';
                case SchedStrategy.LJF
                    text = 'ljf';
                case SchedStrategy.PS
                    text = 'ps';
                case SchedStrategy.DPS
                    text = 'dps';
                case SchedStrategy.GPS
                    text = 'gps';
                case SchedStrategy.SEPT
                    text = 'sept';
                case SchedStrategy.LEPT
                    text = 'lept';
                case SchedStrategy.SRPT
                    text = 'srpt';
                case SchedStrategy.SRPTPRIO
                    text = 'srptprio';
                case SchedStrategy.HOL
                    text = 'hol';
                case SchedStrategy.FCFSPRIO
                    text = 'hol';
                case SchedStrategy.FORK
                    text = 'fork';
                case SchedStrategy.EXT
                    text = 'ext';
                case SchedStrategy.REF
                    text = 'ref';
                case SchedStrategy.LCFSPR
                    text = 'lcfspr';
                case SchedStrategy.LCFSPI
                    text = 'lcfspi';
                case SchedStrategy.LCFSPRIO
                    text = 'lcfsprio';
                case SchedStrategy.LCFSPRPRIO
                    text = 'lcfsprprio';
                case SchedStrategy.LCFSPIPRIO
                    text = 'lcfspiprio';
                case SchedStrategy.FCFSPRPRIO
                    text = 'fcfsprprio';
                case SchedStrategy.FCFSPIPRIO
                    text = 'fcfspiprio';
                case SchedStrategy.POLLING
                    text = 'polling';
                case SchedStrategy.PSPRIO
                    text = 'psprio';
                case SchedStrategy.DPSPRIO
                    text = 'dpsprio';
                case SchedStrategy.GPSPRIO
                    text = 'gpsprio';
                case SchedStrategy.EDD
                    text = 'edd';
                case SchedStrategy.EDF
                    text = 'edf';
                case SchedStrategy.LPS
                    text = 'lps';
                case SchedStrategy.PSJF
                    text = 'psjf';
                case SchedStrategy.FB
                    text = 'fb';
                case SchedStrategy.LRPT
                    text = 'lrpt';
                case SchedStrategy.SETF
                    text = 'setf';
                otherwise
                    line_error(mfilename, 'Unrecognized scheduling strategy type.');
            end
        end

        function type = fromText(text)
            % TYPE = FROMTEXT(TEXT)
            if iscell(text)
                text = text{:};
            end
            switch lower(text)
                case 'inf'
                    type = SchedStrategy.INF;
                case 'fcfs'
                    type = SchedStrategy.FCFS;
                case 'fcfspr'
                    type = SchedStrategy.FCFSPR;
                case 'fcfspi'
                    type = SchedStrategy.FCFSPI;
                case 'lcfs'
                    type = SchedStrategy.LCFS;
                case 'siro'
                    type = SchedStrategy.SIRO;
                case 'sjf'
                    type = SchedStrategy.SJF;
                case 'ljf'
                    type = SchedStrategy.LJF;
                case 'ps'
                    type = SchedStrategy.PS;
                case 'dps'
                    type = SchedStrategy.DPS;
                case 'gps'
                    type = SchedStrategy.GPS;
                case 'sept'
                    type = SchedStrategy.SEPT;
                case 'lept'
                    type = SchedStrategy.LEPT;
                case 'srpt'
                    type = SchedStrategy.SRPT;
                case 'srptprio'
                    type = SchedStrategy.SRPTPRIO;
                case 'hol'
                    type = SchedStrategy.HOL;
                case 'fcfsprio'
                    type = SchedStrategy.FCFSPRIO;
                case 'fork'
                    type = SchedStrategy.FORK;
                case 'ext'
                    type = SchedStrategy.EXT;
                case 'ref'
                    type = SchedStrategy.REF;
                case 'lcfspr'
                    type = SchedStrategy.LCFSPR;
                case 'lcfspi'
                    type = SchedStrategy.LCFSPI;
                case 'lcfsprio'
                    type = SchedStrategy.LCFSPRIO;
                case 'lcfsprprio'
                    type = SchedStrategy.LCFSPRPRIO;
                case 'lcfspiprio'
                    type = SchedStrategy.LCFSPIPRIO;
                case 'fcfsprprio'
                    type = SchedStrategy.FCFSPRPRIO;
                case 'fcfspiprio'
                    type = SchedStrategy.FCFSPIPRIO;
                case 'polling'
                    type = SchedStrategy.POLLING;
                case 'psprio'
                    type = SchedStrategy.PSPRIO;
                case 'dpsprio'
                    type = SchedStrategy.DPSPRIO;
                case 'gpsprio'
                    type = SchedStrategy.GPSPRIO;
                case 'edd'
                    type = SchedStrategy.EDD;
                case 'edf'
                    type = SchedStrategy.EDF;
                case 'lps'
                    type = SchedStrategy.LPS;
                case 'psjf'
                    type = SchedStrategy.PSJF;
                case 'fb'
                    type = SchedStrategy.FB;
                case 'las'
                    type = SchedStrategy.LAS;
                case 'lrpt'
                    type = SchedStrategy.LRPT;
                case 'setf'
                    type = SchedStrategy.SETF;
                case 'set'
                    type = SchedStrategy.SET;
                case 'pp'
                    % LQNS preemptive priority - maps to FCFS with preemptive resume priority
                    type = SchedStrategy.FCFSPRPRIO;
                case 'cfs'
                    % LQNS completely fair scheduling - maps to Generalized Processor Sharing
                    type = SchedStrategy.GPS;
                otherwise
                    line_error(mfilename, 'Unrecognized scheduling strategy text.');
            end
        end

        function property = toProperty(text)
            % PROPERTY = TOPROPERTY(TEXT)            
            if iscell(text)
                text = text{:};            
            end
            switch lower(text)
                case 'inf'
                    property = 'INF';
                case 'fcfs'
                    property = 'FCFS';
                case 'fcfspr'
                    property = 'FCFSPR';
                case 'fcfspi'
                    property = 'FCFSPI';
                case 'lcfs'
                    property = 'LCFS';
                case 'siro'
                    property = 'SIRO';
                case 'sjf'
                    property = 'SJF';
                case 'ljf'
                    property = 'LJF';
                case 'ps'
                    property = 'PS';
                case 'dps'
                    property = 'DPS';
                case 'gps'
                    property = 'GPS';
                case 'sept'
                    property = 'SEPT';
                case 'lept'
                    property = 'LEPT';
                case 'srpt'
                    property = 'SRPT';
                case 'srptprio'
                    property = 'SRPTPRIO';
                case 'hol'
                    property = 'HOL';
                case 'fcfsprio'
                    property = 'FCFSPRIO';
                case 'fork'
                    property = 'FORK';
                case 'ext'
                    property = 'EXT';
                case 'ref'
                    property = 'REF';
                case 'lcfspr'
                    property = 'LCFSPR';
                case 'lcfspi'
                    property = 'LCFSPI';
                case 'lcfsprio'
                    property = 'LCFSPRIO';
                case 'lcfsprprio'
                    property = 'LCFSPRPRIO';
                case 'lcfspiprio'
                    property = 'LCFSPIPRIO';
                case 'fcfsprprio'
                    property = 'FCFSPRPRIO';
                case 'fcfspiprio'
                    property = 'FCFSPIPRIO';
                case 'polling'
                    property = 'POLLING';
                case 'psprio'
                    property = 'PSPRIO';
                case 'dpsprio'
                    property = 'DPSPRIO';
                case 'gpsprio'
                    property = 'GPSPRIO';
                case 'edd'
                    property = 'EDD';
                case 'edf'
                    property = 'EDF';
                case 'lps'
                    property = 'LPS';
                case 'psjf'
                    property = 'PSJF';
                case 'fb'
                    property = 'FB';
                case 'las'
                    property = 'LAS';
                case 'lrpt'
                    property = 'LRPT';
                case 'setf'
                    property = 'SETF';
                case 'set'
                    property = 'SET';
                otherwise
                    line_error(mfilename, 'Unrecognized scheduling strategy property.');
            end
        end

        function text = toFeature(type)
            % TEXT = TOFEATURE(TYPE)
            switch type
                case SchedStrategy.INF
                    text = 'SchedStrategy_INF';
                case SchedStrategy.FCFS
                    text = 'SchedStrategy_FCFS';
                case SchedStrategy.FCFSPR
                    text = 'SchedStrategy_FCFSPR';
                case SchedStrategy.FCFSPI
                    text = 'SchedStrategy_FCFSPI';
                case SchedStrategy.LCFS
                    text = 'SchedStrategy_LCFS';
                case SchedStrategy.SIRO
                    text = 'SchedStrategy_SIRO';
                case SchedStrategy.SJF
                    text = 'SchedStrategy_SJF';
                case SchedStrategy.LJF
                    text = 'SchedStrategy_LJF';
                case SchedStrategy.PS
                    text = 'SchedStrategy_PS';
                case SchedStrategy.DPS
                    text = 'SchedStrategy_DPS';
                case SchedStrategy.GPS
                    text = 'SchedStrategy_GPS';
                case SchedStrategy.SEPT
                    text = 'SchedStrategy_SEPT';
                case SchedStrategy.LEPT
                    text = 'SchedStrategy_LEPT';
                case SchedStrategy.SRPT
                    text = 'SchedStrategy_SRPT';
                case SchedStrategy.SRPTPRIO
                    text = 'SchedStrategy_SRPTPRIO';
                case SchedStrategy.HOL
                    text = 'SchedStrategy_HOL';
                case SchedStrategy.FCFSPRIO
                    text = 'SchedStrategy_HOL';
                case SchedStrategy.FORK
                    text = 'SchedStrategy_FORK';
                case SchedStrategy.EXT
                    text = 'SchedStrategy_EXT';
                case SchedStrategy.REF
                    text = 'SchedStrategy_REF';
                case SchedStrategy.LCFSPR
                    text = 'SchedStrategy_LCFSPR';
                case SchedStrategy.LCFSPI
                    text = 'SchedStrategy_LCFSPI';
                case SchedStrategy.LCFSPRIO
                    text = 'SchedStrategy_LCFSPRIO';
                case SchedStrategy.LCFSPRPRIO
                    text = 'SchedStrategy_LCFSPRPRIO';
                case SchedStrategy.LCFSPIPRIO
                    text = 'SchedStrategy_LCFSPIPRIO';
                case SchedStrategy.FCFSPRPRIO
                    text = 'SchedStrategy_FCFSPRPRIO';
                case SchedStrategy.FCFSPIPRIO
                    text = 'SchedStrategy_FCFSPIPRIO';
                case SchedStrategy.POLLING
                    text = 'SchedStrategy_POLLING';
                case SchedStrategy.PSPRIO
                    text = 'SchedStrategy_PSPRIO';
                case SchedStrategy.DPSPRIO
                    text = 'SchedStrategy_DPSPRIO';
                case SchedStrategy.GPSPRIO
                    text = 'SchedStrategy_GPSPRIO';
                case SchedStrategy.EDD
                    text = 'SchedStrategy_EDD';
                case SchedStrategy.EDF
                    text = 'SchedStrategy_EDF';
                case SchedStrategy.LPS
                    text = 'SchedStrategy_LPS';
                case SchedStrategy.PSJF
                    text = 'SchedStrategy_PSJF';
                case SchedStrategy.FB
                    text = 'SchedStrategy_FB';
                case SchedStrategy.LRPT
                    text = 'SchedStrategy_LRPT';
                case SchedStrategy.SETF
                    text = 'SchedStrategy_SETF';
                otherwise
                    line_error(mfilename, 'Unrecognized scheduling strategy feature.');
            end
        end
    end
end
