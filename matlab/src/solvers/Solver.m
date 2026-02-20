classdef Solver < handle
    % Abstract base class for all LINE model solution algorithms
    %
    % Provides common interface and infrastructure for queueing network analysis.
    %
    % Copyright (c) 2012-2026, Imperial College London
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
            self.model = model; % passed by reference for cache update
            self.name = name;
            self.options = options;
            self.enableChecks = true;
        end
    end

    methods (Abstract) % implemented with errors for Octave compatibility
        sn = getStruct(self);
        runtime = runAnalyzer(self, options)% generic method to run the solver
        bool = supports(self,model);
    end


    methods
        function model = getModel(self)
            model = self.model;
        end

        function self = setChecks(self, bool)
            self.enableChecks = bool;
        end
        function out = getName(self)
            % OUT = GETNAME()
            % Get solver name
            out = self.name;
        end
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
            self.result = [];
            Solver.resetRandomGeneratorSeed(self.options.seed);
        end

        function self = setOptions(self, options)
            % SELF = SETOPTIONS(OPTIONS)
            % Set a new options data structure
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

        function [optList, allOpt] = listValidOptions()
            % OPTLIST = LISTVALIDOPTIONS()
            % List valid fields for options data structure
            optList = {'cache','cutoff','force','init_sol','iter_max','iter_tol','lang','tol', ...
                'keep','method','odesolvers','samples','seed','stiff', 'timespan','timestep','verbose','config.multiserver','config.fork_join','fork_join','confint'};
            %,'amva.aql','aql','amva.qdaql','qdaql','mom','brute','jmva.ls','jmt.jmva.ls'
            allOpt = {'cache','cutoff','force','init_sol','iter_max','iter_tol','lang','tol', ...
                'keep','method','odesolvers','samples','seed','stiff', 'timespan','timestep','verbose','config.multiserver','config.fork_join','fork_join','confint', ...
                'default','exact','auto','ctmc','ctmc.gpu','gpu','mva','mva.exact','mva.amva','mva.qna','sqni','mva.sqni',...
                'amva','amva.bs','amva.qd','bs','qd','amva.qli','qli','amva.fli','fli','amva.lin','lin','amva.qdlin','qdlin',...
                'ssa','ssa.parallel','serial','parallel','nrm',...
                'jmt','jsim','replication','jmva','jmva.amva','jmva.mva','jmva.recal','jmva.mom','jmva.comom','jmva.chow','jmva.bs','jmva.aql','jmva.lin','jmva.dmlin','jmt.jsim',...
                'jmt.jmva','jmt.jmva.mva','jmt.jmva.amva','jmt.jmva.recal','jmt.jmva.comom','jmt.jmva.chow','jmt.jmva.bs','jmt.jmva.aql','jmt.jmva.lin','jmt.jmva.dmlin',...
                'ca','comom','comomld','gm','propfair','recal','kt', 'rd', 'nrp', 'nrl', ...
                'nc.brute','nc.ca','nc.comom','nc.comomld','nc.gm','nc.mom','nc.propfair','nc.recal','nc.kt', 'nc.rd', 'nc.nr.probit', 'nc.nr.logit', ...
                'fluid','matrix','softmin','statedep','closing','diffusion','fluid.softmin','fluid.statedep','fluid.closing','fluid.matrix','fluid.diffusion',...
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
            options = lineDefaults;
        end

        function [enabled, level] = parseConfInt(confint)
            % [ENABLED, LEVEL] = PARSECONFINT(CONFINT)
            % Parse confidence interval option value
            % Returns enabled (true/false) and confidence level (0.0-1.0)
            if islogical(confint)
                enabled = confint;
                level = 0.95;  % default 95% confidence
            elseif isnumeric(confint) && confint > 0 && confint < 1
                enabled = true;
                level = confint;
            elseif isnumeric(confint) && confint == 0
                enabled = false;
                level = 0.95;
            else
                enabled = false;
                level = 0.95;
            end
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
                elseif isscalar(varargin)
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
                                    case 'fork_join'
                                        options.config.fork_join = varargin{2};
                                    case 'verbose'
                                        % Convert string verbose level to VerboseLevel constant
                                        if ischar(varargin{2}) || isstring(varargin{2})
                                            switch upper(char(varargin{2}))
                                                case 'SILENT'
                                                    options.verbose = VerboseLevel.SILENT;
                                                case 'STD'
                                                    options.verbose = VerboseLevel.STD;
                                                case 'DEBUG'
                                                    options.verbose = VerboseLevel.DEBUG;
                                                otherwise
                                                    options.verbose = varargin{2};
                                            end
                                        else
                                            options.verbose = varargin{2};
                                        end
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

    end
end
