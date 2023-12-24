classdef (Sealed) SchedStrategy
    % Enumeration of scheduling strategies
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        INF = 'inf'; % infinite server
        FCFS = 'fcfs';
        LCFS = 'lcfs';
        LCFSPR = 'lcfspr';
        SIRO = 'siro'; % service in random order
        SJF = 'sjf';
        LJF = 'ljf';
        PS = 'ps'; % egalitarian PS
        DPS = 'dps';
        GPS = 'gps';
        SEPT = 'sept';
        LEPT = 'lept';
        HOL = 'hol';
        FORK = 'fork';
        EXT = 'ext'; % external world (open arrival source and sink)
        REF = 'ref'; % reference node in LayeredNetworks
        POLLING = 'polling';
        
        ID_INF = 0;
        ID_FCFS = 1;
        ID_LCFS = 2;
        ID_SIRO = 3;
        ID_SJF = 4;
        ID_LJF = 5;
        ID_PS = 6;
        ID_DPS = 7;
        ID_GPS = 8;
        ID_SEPT = 9;
        ID_LEPT = 10;
        ID_HOL = 11;
        ID_FORK = 12;
        ID_EXT = 13;
        ID_REF = 14;
        ID_LCFSPR = 15;
        ID_POLLING = 16;
    end
    
    methods (Static)
        function id = toId(type)
            % ID = TOID(TYPE)
            switch type
                case SchedStrategy.INF
                    id = 0;
                case SchedStrategy.FCFS
                    id = 1;
                case SchedStrategy.LCFS
                    id = 2;
                case SchedStrategy.SIRO
                    id = 3;
                case SchedStrategy.SJF
                    id = 4;
                case SchedStrategy.LJF
                    id = 5;
                case SchedStrategy.PS
                    id = 6;
                case SchedStrategy.DPS
                    id = 7;
                case SchedStrategy.GPS
                    id = 8;
                case SchedStrategy.SEPT
                    id = 9;
                case SchedStrategy.LEPT
                    id = 10;
                case SchedStrategy.HOL
                    id = 11;
                case SchedStrategy.FORK
                    id = 12;
                case SchedStrategy.EXT
                    id = 13;
                case SchedStrategy.REF
                    id = 14;
                case SchedStrategy.LCFSPR
                    id = 15;
                case SchedStrategy.POLLING
                    id = 16;
            end
        end

  function type = fromId(id)
            % TYPE = FROMID(ID)
            switch id
                case SchedStrategy.ID_INF
                    type=SchedStrategy.INF;
                case SchedStrategy.ID_FCFS
                    type=SchedStrategy.FCFS;
                case SchedStrategy.ID_LCFS
                    type=SchedStrategy.LCFS;
                case SchedStrategy.ID_SIRO
                    type=SchedStrategy.SIRO;
                case SchedStrategy.ID_SJF
                    type=SchedStrategy.SJF;
                case SchedStrategy.ID_LJF
                    type=SchedStrategy.LJF;
                case SchedStrategy.ID_PS
                    type=SchedStrategy.PS;
                case SchedStrategy.ID_DPS
                    type=SchedStrategy.DPS;
                case SchedStrategy.ID_GPS
                    type=SchedStrategy.GPS;
                case SchedStrategy.ID_SEPT
                    type=SchedStrategy.SEPT;
                case SchedStrategy.ID_LEPT
                    type=SchedStrategy.LEPT;
                case SchedStrategy.ID_HOL
                    type=SchedStrategy.HOL;
                case SchedStrategy.ID_FORK
                    type=SchedStrategy.FORK;
                case SchedStrategy.ID_EXT
                    type=SchedStrategy.EXT;
                case SchedStrategy.ID_REF
                    type=SchedStrategy.REF;
                case SchedStrategy.ID_LCFSPR
                    type=SchedStrategy.LCFSPR;
                case SchedStrategy.ID_POLLING
                    type=SchedStrategy.POLLING;
            end
        end
        
        function text = toText(type)
            % TEXT = TOID(TYPE)
            switch type
                case SchedStrategy.INF
                    text = 'inf';
                case SchedStrategy.FCFS
                    text = 'fcfs';
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
                case SchedStrategy.HOL
                    text = 'hol';
                case SchedStrategy.FORK
                    text = 'fork';
                case SchedStrategy.EXT
                    text = 'ext';
                case SchedStrategy.REF
                    text = 'ref';
                case SchedStrategy.LCFSPR
                    text = 'lcfspr';
                case SchedStrategy.POLLING
                    text = 'polling';
            end
        end
        
   function type = fromText(text)
            % TEXT = TOID(TYPE)
            switch text
                case 'inf'
                    type = SchedStrategy.INF;
                case 'fcfs'
                    type = SchedStrategy.FCFS;
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
                case 'hol'
                    type = SchedStrategy.HOL;
                case 'fork'
                    type = SchedStrategy.FORK;
                case 'ext'
                    type = SchedStrategy.EXT;
                case 'ref'
                    type = SchedStrategy.REF;
                case 'lcfspr'
                    type = SchedStrategy.LCFSPR;
                case 'polling'
                    type = SchedStrategy.POLLING;
            end
        end        
        
        function property = toProperty(text)
            % PROPERTY = TOPROPERTY(TEXT)            
            if iscell(text)
                text=text{:};
            end
            switch text
                case 'inf'
                    property = 'INF';
                case 'fcfs'
                    property = 'FCFS';
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
                case 'hol'
                    property = 'HOL';
                case 'ext'
                    property = 'EXT';
                case 'fork'
                    property = 'FORK';
                case 'ref'
                    property = 'REF';
                case 'ext'
                    property = 'EXT';
                case 'lcfspr'
                    property = 'LCFSPR';
                case 'polling'
                    property = 'POLLING';
            end
        end
        
        function text = toFeature(type)
            % TEXT = TOFEATURE(TYPE)
            
            switch type
                case SchedStrategy.INF
                    text = 'SchedStrategy_INF';
                case SchedStrategy.FCFS
                    text = 'SchedStrategy_FCFS';
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
                case SchedStrategy.HOL
                    text = 'SchedStrategy_HOL';
                case SchedStrategy.FORK
                    text = 'SchedStrategy_FORK';
                case SchedStrategy.EXT
                    text = 'SchedStrategy_EXT';
                case SchedStrategy.REF
                    text = 'SchedStrategy_REF';
                case SchedStrategy.LCFSPR
                    text = 'SchedStrategy_LCFSPR';
                case SchedStrategy.POLLING
                    text = 'SchedStrategy_POLLING';
            end
        end
    end
    
end
