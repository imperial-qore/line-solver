classdef (Sealed) ProcessType
    % Enumeration of process ts
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties (Constant)       
        EXP = 0;
        ERLANG = 1;
        HYPEREXP = 2;
        PH = 3;
        APH = 4;
        MAP = 5;
        UNIFORM = 6;
        DET = 7;
        COXIAN = 8;
        GAMMA = 9;
        PARETO = 10;
        MMPP2 = 11;
        REPLAYER = 12;
        TRACE = 12;
        IMMEDIATE = 13;
        DISABLED = 14;
        COX2 = 15;
        WEIBULL = 16;
        LOGNORMAL = 17;
        DUNIFORM = 18;
        BERNOULLI = 19;
        PRIOR = 20;
        BINOMIAL = 21;
        POISSON = 22;
        GEOMETRIC = 23;
        BMAP = 24;
        ME = 25;
        RAP = 26;
        DISCRETESAMPLER = 27;
        ZIPF = 28;
        DMAP = 29;
        MMAP = 31;
        EMPIRICALCDF = 32;
        NHPP = 33;
    end
    
    methods (Static)
        function bool = isMarkovian(t)
            % BOOL = ISMARKOVIAN(T)
            %
            % True when sn.proc carries an exact matrix representation of a
            % process of type T: a genuine (D0,D1) pair, or its
            % matrix-exponential analogue for ME/RAP.
            %
            % The distinction is about SN.PROC, not about what getProcess
            % returns. getProcess hands back raw distribution PARAMETERS for
            % several non-Markovian families -- Gamma, Weibull, Lognormal,
            % Pareto and Uniform return two scalars (Pareto {alpha,k}, Uniform
            % {min,max}) -- and refreshProcessRepresentations replaces those
            % with convertToMAP = map_erlang(mean, n) before storing them in
            % sn.proc, where n = ceil(1/SCV) capped at 100 (n = 20 when
            % SCV < CoarseTol). That fit matches the mean, and matches the SCV
            % only when SCV <= 1: Pareto with SCV 64 gives n = 1, a single
            % exponential of SCV 1. So for these types sn.proc is an
            % approximation, not the law that was requested, and nothing on the
            % cell says so -- the only signal is sn.procid.
            %
            % Solvers that read sn.proc as if it were the exact law must gate on
            % this predicate. It is the procid-level counterpart of the JAR's
            % Distribution.isMarkovian(), i.e. of the Markovian class hierarchy,
            % so the two lists must stay in step.
            bool = any(t == [ProcessType.EXP, ProcessType.ERLANG, ...
                ProcessType.HYPEREXP, ProcessType.PH, ProcessType.APH, ...
                ProcessType.MAP, ProcessType.COXIAN, ProcessType.COX2, ...
                ProcessType.MMPP2, ProcessType.ME, ProcessType.RAP, ...
                ProcessType.DMAP, ProcessType.BMAP, ProcessType.MMAP]);
        end

        function t = fromId(id)
            % ID = TOID(TYPE)
            t = id;
        end
        
        function id = toId(t)
            % ID = TOID(TYPE)            
            id = t;
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
                case {'Replayer', 'Trace'}
                    t = ProcessType.REPLAYER;
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
                case 'DiscreteUniform'
                    t = ProcessType.DUNIFORM;
                case 'Bernoulli'
                    t = ProcessType.BERNOULLI;
                case 'Prior'
                    t = ProcessType.PRIOR;
                case 'Binomial'
                    t = ProcessType.BINOMIAL;
                case 'Poisson'
                    t = ProcessType.POISSON;
                case 'Geometric'
                    t = ProcessType.GEOMETRIC;
                case 'BMAP'
                    t = ProcessType.BMAP;
                case {'ME', 'CME'}
                    % A CME is an ME: the concentrated representation is a
                    % subclass of ME and carries no distinct process type, so
                    % every solver gate and sn.procid entry that accepts ME
                    % accepts it. refreshProcessTypes passes the class name, so
                    % the subclass has to be listed here explicitly.
                    t = ProcessType.ME;
                case 'RAP'
                    t = ProcessType.RAP;
                case 'DiscreteSampler'
                    t = ProcessType.DISCRETESAMPLER;
                case 'Zipf'
                    t = ProcessType.ZIPF;
                case 'DMAP'
                    t = ProcessType.DMAP;
                case 'NHPP'
                    t = ProcessType.NHPP;
                case 'EmpiricalCDF'
                    % Class name, as passed by refreshProcessTypes. Without this
                    % case an EmpiricalCDF service aborts fromText during
                    % getStruct, before any solver feature gate can report that
                    % the distribution is unsupported.
                    t = ProcessType.EMPIRICALCDF;
                case {'MarkedMAP', 'MMAP', 'MarkedMMPP'}
                    t = ProcessType.MMAP;
                otherwise
                    line_error(mfilename, sprintf('Unrecognized process type: %s', text));
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
                case ProcessType.DUNIFORM
                    text = 'DiscreteUniform';
                case ProcessType.BERNOULLI
                    text = 'Bernoulli';
                case ProcessType.PRIOR
                    text = 'Prior';
                case ProcessType.BINOMIAL
                    text = 'Binomial';
                case ProcessType.POISSON
                    text = 'Poisson';
                case ProcessType.GEOMETRIC
                    text = 'Geometric';
                case ProcessType.BMAP
                    text = 'BMAP';
                case ProcessType.ME
                    text = 'ME';
                case ProcessType.RAP
                    text = 'RAP';
                case ProcessType.DISCRETESAMPLER
                    text = 'DiscreteSampler';
                case ProcessType.ZIPF
                    text = 'Zipf';
                case ProcessType.DMAP
                    text = 'DMAP';
                case ProcessType.NHPP
                    text = 'NHPP';
                case ProcessType.MMAP
                    text = 'MMAP';
                case ProcessType.EMPIRICALCDF
                    text = 'EmpiricalCDF';
                otherwise
                    if isnan(t)
                        text = 'Unknown';
                    else
                        text = sprintf('Unknown(%d)', t);
                    end
            end

        end

    end
end
