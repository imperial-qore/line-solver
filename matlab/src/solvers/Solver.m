classdef Solver < handle
    % Abstract class for model solution algorithms and tools
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.

    properties (Constant)
        VERBOSE_SILENT = 0;
        VERBOSE_STD = 1;
        VERBOSE_DEBUG = 2;
    end

    properties (Hidden, Access = public)
        enableChecks;
    end

    properties (Access = public)
        options; % Data structure with solver options
        name; % Solver name
        model; % Model to be solved
        result; % last result
        obj;
    end

    methods (Hidden)
        %Constructor
        function self = Solver(model, name, options)
            % SELF = SOLVER(MODEL, NAME, OPTIONS)
            if nargin<3 %~exist('options','var')
                options = self.defaultOptions;
            end
            self.model = model;
            self.name = name;
            self.options = options;
            self.enableChecks = true;
        end
    end

    methods (Static)
        function bool = supports(self,model)
            % BOOL = SUPPORTS(SELF,MODEL)
            % True if the input model is supported by the solver
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end
    end

    methods %(Abstract) % implemented with errors for Octave compatibility

        function model = getModel(self)
            model = self.model;
        end

        sn = getStruct(self);

        function runtime = runAnalyzer(self, options) % generic method to run the solver
            % RUNTIME = RUNANALYZER()
            % Run the solver % GENERIC METHOD TO RUN THE SOLVER
            % Solve the model
            line_error(mfilename,'An abstract method was called. The function needs to be overridden by a subclass.');
        end

        function self = setChecks(self, bool)
            self.enableChecks = bool;
        end
    end

    methods
        function out = getName(self)
            % OUT = GETNAME()
            % Get solver name
            out = self.name;
        end

        function options = getDefaultOptions(self)
            % OPTIONS = GETDEFAULTOPTIONS()
            % Get option data structure with default values
            options = Solver.defaultOptions;
        end

    end

    methods

        function results = getResults(self)
            % RESULTS = GETRESULTS()
            % Return results data structure
            results = self.result;
        end

        function bool = hasResults(self)
            % BOOL = HASRESULTS()
            % Check if the model has been solved
            bool = ~isempty(self.result);
        end

        function options = getOptions(self)
            % OPTIONS = GETOPTIONS()
            % Return options data structure
            options = self.options;
        end

        function reset(self)
            % RESET()
            % Dispose previously stored results
            resetResults(self);
        end

        function resetResults(self)
            % RESET()
            % Dispose previously stored results
            self.result = [];
        end

        function self = setOptions(self, options)
            % SELF = SETOPTIONS(OPTIONS)
            % Set a new options data structure
            defaultOptions = self.defaultOptions;
            optList = Solver.listValidOptions();
            for l=1:length(optList)
                if ~isfield(options,optList{l})
                    switch optList{l}
                        case 'config.multiserver'
                            % no-op
                        otherwise
                            options.(optList{l}) = defaultOptions.(optList{l});
                    end
                end
            end
            self.options = options;
        end

    end

    methods (Static)

        function resetRandomGeneratorSeed(seed)
            % RESETRANDOMGENERATORSEED(SEED)
            % Assign a new seed to the random number generator
            warning('off','MATLAB:RandStream:ActivatingLegacyGenerators');
            warning('off','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState');
            rng(seed,'twister');
        end

        function bool = isAvailable()
            % BOOL = ISAVAILABLE()
            % Check if external dependencies are available for the solver
            bool = true;
        end

        function bool = isJavaAvailable()
            % BOOL = ISJAVAAVAILABLE()
            % Check if Java dependencies are available for the solver
            bool = true;
            if ispc % windows
                [~,ret] = dos('java -version');
                if strfind(ret,'not recognized') %#ok<STRIFCND>
                    bool = false;
                end
            else %linux
                [~,ret] = unix('java -version');
                if strfind(ret,'command not found') %#ok<STRIFCND>
                    bool = false;
                end
            end
        end

        function fun = accurateStiffOdeSolver()
            % FUN = ACCURATESTIFFODESOLVER()
            % Return default high-accuracy stiff solver
            fun = @ode15s; % variable-order method
        end

        function fun = accurateOdeSolver()
            % FUN = ACCURATEODESOLVER()
            % Return default high-accuracy non-stiff solver
            fun = @ode113; % variable-order method
        end

        function fun = fastStiffOdeSolver()
            % FUN = FASTSTIFFODESOLVER()
            % Return default low-accuracy stiff solver
            fun = @ode23s; % low-order method
        end

        function fun = fastOdeSolver()
            % FUN = FASTODESOLVER()
            % Return default low-accuracy non-stiff solver
            fun = @ode23; % low-order method
        end

        function [optList, allOpt] = listValidOptions()
            % OPTLIST = LISTVALIDOPTIONS()
            % List valid fields for options data structure
            optList = {'cache','cutoff','force','init_sol','iter_max','iter_tol','lang','tol', ...
                'keep','method','odesolvers','samples','seed','stiff', 'timespan','verbose','config.multiserver'};
            allOpt = {'cache','cutoff','force','init_sol','iter_max','iter_tol','lang','tol', ...
                'keep','method','odesolvers','samples','seed','stiff', 'timespan','verbose','config.multiserver', ...
                'default','exact','auto','ctmc','ctmc.gpu','gpu','mva','mva.exact','mva.amva','mva.qna','sqrt','mva.sqrt',...
                'amva','amva.bs','amva.qd','bs','qd','amva.qli','qli','amva.fli','fli','amva.aql','aql','amva.qdaql','qdaql','amva.lin','lin','amva.qdlin','qdlin',...
                'ssa','ssa.serial.hash','ssa.para.hash','ssa.parallel.hash','ssa.serial','ssa.para','ssa','taussa','tauleap',...
                'ssa.parallel','serial.hash','serial','para','parallel','para.hash','parallel.hash',...
                'jmt','jsim','jmva','jmva.amva','jmva.mva','jmva.recal','jmva.mom','jmva.comom','jmva.chow','jmva.bs','jmva.aql','jmva.lin','jmva.dmlin','jmva.ls',...
                'jmt.jsim','jmt.jmva','jmt.jmva.mva','jmt.jmva.amva','jmt.jmva.recal','jmt.jmva.comom','jmt.jmva.chow','jmt.jmva.bs','jmt.jmva.aql','jmt.jmva.lin','jmt.jmva.dmlin','jmt.jmva.ls',...
                'brute','ca','comom','comomld','gm','mom','propfair','recal','kt', 'rd', 'nr.probit', 'nr.logit', 'nc.brute','nc.ca','nc.comom','nc.comomld','nc.gm','nc.mom','nc.propfair','nc.recal','nc.kt', 'nc.rd', 'nc.nr.probit', 'nc.nr.logit', ...
                'fluid','matrix','softmin','statedep','closing','fluid.softmin','fluid.statedep''fluid.closing','fluid.matrix',...
                'nc','nc.exact','nc.imci','ls','nc.ls','nc.cub','cub','le','nc.le','nc.panacea','panacea','nc.mmint2','mmint2','nc.gleint','gleint','mam','dec.source','dec.mmap',...
                'mmk','gigk', 'gigk.kingman_approx', ...
                'mm1','mg1','gm1','gig1','gim1','gig1.kingman','gig1.gelenbe','gig1.heyman','gig1.kimura','gig1.allen','gig1.kobayashi','gig1.klb','gig1.marchal',...
                'aba.upper','aba.lower','gb.upper','gb.lower','sb.upper','sb.lower','bjb.upper','bjb.lower','pb.upper','pb.lower'};
        end

        function bool = isValidOption(optName)
            % BOOL = ISVALIDOPTION(OPTNAME)
            % Check if the given option exists for the solver
            [~,allOpts] = Solver.listValidOptions();
            bool = any(cell2mat(findstring(optName, allOpts))==1);
        end

        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            % Return default options
            options = lineDefaults('');
            options.iter_max = 100;       
        end

        function options = parseOptions(varargin, defaultOptions)
            % OPTIONS = PARSEOPTIONS(VARARGIN, DEFAULTOPTIONS)
            % Parse option parameters into options data structure
            if isempty(varargin)
                options = defaultOptions;
            elseif isstruct(varargin{1})
                options = varargin{1};
            elseif ischar(varargin{1})
                if length(varargin)>1 && isstruct(varargin{2}) % options struct after method field
                    options = varargin{2};
                    varargin(2) = [];
                elseif length(varargin)==1
                    options = defaultOptions;
                    options.method = varargin{1};
                else
                    options = defaultOptions;
                end
                [optList, allOpt] = Solver.listValidOptions();
                allMethodsList = setdiff(allOpt, optList);
                while ~isempty(varargin)
                    if Solver.isValidOption(varargin{1}) || startsWith(varargin{1},'config.')
                        switch varargin{1}
                            case allMethodsList
                                options.method = varargin{1};
                                varargin(1) = [];
                            otherwise
                                switch varargin{1}
                                    case 'config.eventcache'
                                        options.config.eventcache = varargin{2};
                                    case 'config.highvar'
                                        options.config.highvar = varargin{2};
                                    case 'config.multiserver'
                                        options.config.multiserver = varargin{2};
                                    case 'config.np_priority'
                                        options.config.np_priority = varargin{2};
                                    case 'config.fork_join'
                                        options.config.fork_join = varargin{2};
                                    otherwise
                                        options.(varargin{1}) = varargin{2};
                                end
                                varargin(1) = [];
                                varargin(1) = [];
                        end
                    else
                        %line_warning(mfilename,sprintf('Option "%s" does not exist. Ignoring.',varargin{1}));
                        varargin(1) = [];
                    end
                end
            else
                line_error(mfilename,'Invalid parameter.');
            end
        end

        function solver = load(chosenmethod, model, varargin)
            % SOLVER = LOAD(CHOSENMETHOD, MODEL, VARARGIN)
            % Returns a solver configured to run the chosen method
            options = Solver.parseOptions(varargin, Solver.defaultOptions);
            options.method = chosenmethod;
            switch options.method
                case {'default','auto'}
                    if strcmp(options.method,'auto'), options.method='default'; end
                    solver = LINE(model, options);
                case {'ctmc','ctmc.gpu','gpu'}
                    if strcmp(options.method,'ctmc'), options.method='default'; end
                    options.method = erase(options.method,'ctmc.');
                    solver = SolverCTMC(model, options);
                case {'mva','mva.exact','amva','mva.amva','qna','mva.qna','sqrt','mva.sqrt' ...
                        'amva.qd','mva.amva.qd', ...
                        'amva.bs','mva.amva.bs', ...
                        'amva.qli','mva.amva.qli', ...
                        'amva.fli','mva.amva.fli', ...
                        'amva.aql','mva.amva.aql', ...
                        'amva.lin','mva.amva.lin', ...
                        'amva.qdaql','mva.amva.qdaql', ...
                        'mm1','mmk','mg1','gm1','gig1','gim1','gig1.kingman', ...
                        'gigk', 'gigk.kingman_approx', ...
                        'gig1.gelenbe','gig1.heyman','gig1.kimura','gig1.allen','gig1.kobayashi','gig1.klb','gig1.marchal', ...
                        'aba.upper', 'aba.lower', 'bjb.upper', 'bjb.lower', 'gb.upper', 'gb.lower', 'pb.upper', 'pb.lower', 'sb.upper', 'sb.lower'}
                    if strcmp(options.method,'mva'), options.method='default'; end
                    options.method = erase(options.method,'mva.');
                    solver = SolverMVA(model, options);
                case {'ssa','ssa.serial.hash','ssa.para.hash','ssa.parallel.hash','ssa.serial','ssa.para','ssa.parallel','serial.hash','serial','para','parallel','para.hash','parallel.hash','hashed','taussa','tauleap'}
                    if strcmp(options.method,'ssa'), options.method='default'; end
                    options.method = erase(options.method,'ssa.');
                    solver = SolverSSA(model, options);
                case {'jmt','jsim','jmva','jmva.mva','jmva.recal','jmva.comom','jmva.chow','jmva.bs','jmva.aql','jmva.lin','jmva.dmlin','jmva.ls',...
                        'jmt.jsim','jmt.jmva','jmt.jmva.mva','jmt.jmva.amva','jmva.amva','jmt.jmva.recal','jmt.jmva.comom','jmt.jmva.chow','jmt.jmva.bs','jmt.jmva.aql','jmt.jmva.lin','jmt.jmva.dmlin','jmt.jmva.ls'}
                    if strcmp(options.method,'jmt'), options.method='default'; end
                    options.method = erase(options.method,'jmt.');
                    solver = SolverJMT(model, options);
                case {'fluid','fluid.softmin','fluid.statedep','fluid.closing'}
                    if strcmp(options.method,'fluid'), options.method='default'; end
                    options.method = erase(options.method,'fluid.');
                    solver = SolverFluid(model, options);
                case {'nc','nc.exact','nc.imci','nc.ls','comom','comomld','cub','ls','nc.le','le','mmint2','nc.panacea','nc.pana','nc.mmint2','nc.kt','nc.deterministic','nc.sampling','nc.propfair','nc.comom','nc.comomld','nc.mom','nc.cub','nc.brute','nc.rd', 'nc.nr.probit', 'nc.nr.logit','nc.gm'}
                    if strcmp(options.method,'nc'), options.method='default'; end
                    options.method = erase(options.method,'nc.');
                    solver = SolverNC(model, options);
                case {'mam','mam.dec.source','mam.dec.mmap','mam.dec.poisson'}
                    if strcmp(options.method,'mam'), options.method='default'; end
                    options.method = erase(options.method,'mam.');
                    solver = SolverMAM(model, options);
                otherwise
                    if strcmp(options.method,'auto'), options.method='default'; end
                    solver = LINE(model, options);
            end
        end

    end
end
