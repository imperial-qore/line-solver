classdef SolverMAM < NetworkSolver
    % Matrix-Analytic and RCAT methods solver
    %
    % Implements matrix-analytic methods and RCAT for structured Markov chain analysis.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods
        function self = SolverMAM(model,varargin)
            % SOLVERMAM Create a Matrix-Analytic Methods solver instance
            %
            % @brief Creates a MAM solver for structured Markov chain analysis
            % @param model Network model to be analyzed via matrix-analytic methods
            % @param varargin Optional parameters (method, tolerance, etc.)
            % @return self SolverMAM instance configured for MAM analysis
            
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, self.defaultOptions));
            self.setLang();
        end
                
        function sn = getStruct(self)
            % QN = GETSTRUCT()
            
            % Get data structure summarizing the model
            sn = self.model.getStruct(true);
        end
        
        runtime = runAnalyzer(self, options);
        RD = getCdfRespT(self, R);

            function [allMethods] = listValidMethods(self)
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            sn = self.model.getStruct();
            % Note: Method order must match expected values in test files.
            % New methods should be added at the end to preserve index alignment.
            % 'exact' method removed - autocat moved to line-legacy.git
            allMethods = {'default','dec.source','dec.mmap','dec.poisson','mna','inap','ldqbd','inapplus','dec.source.mmap','inapinf'};
        end

        function featSupported = getMethodFeatureSet(self, method) %#ok<INUSD>
            % Every MAM method shares the solver-level feature envelope; the
            % genuine per-method restrictions ('mna', 'ldqbd') are structural
            % and are applied in supportsModelMethod below.
            %
            % Defining this is what lets NetworkSolver.supportsModelMethod name
            % the offending features: with no method feature set it falls back
            % to the coarse supports(model) and returns an empty reason, so the
            % gate could only report "features not supported" without saying
            % which ones.
            featSupported = SolverMAM.getFeatureSet();
        end

        function [bool, reason] = supportsModelMethod(self, method)
            % Method-aware gate for the genuine per-method restrictions of the
            % MAM analyzer (mirrors the inline guards in solver_mam_analyzer):
            % 'mna' does not support mixed open/closed models, and 'ldqbd'
            % requires a single-class model. All other methods rely on the
            % coarse MAM feature set and the analyzer's topology-based routing
            % (e.g. a Fork-Join model on 'default'/'dec.source' is routed to the
            % FJ solver, not rejected).
            sn = self.model.getStruct();
            switch method
                case 'mna'
                    if ~sn_is_open_model(sn) && ~sn_is_closed_model(sn)
                        bool = false;
                        reason = 'The mna method does not support mixed open/closed models.';
                        return;
                    end
                    % mna rejects self-looping classes (no inter-station flow);
                    % see _kb/06-solver-catalog.md for rationale
                    Vmna = cellsum(sn.visits);
                    for kmna = 1:sn.nclasses
                        % Only a CLOSED class is self-looping (njobs<Inf); see _kb/06-solver-catalog.md for rationale
                        if ~isfinite(sn.njobs(kmna))
                            continue;
                        end
                        vis = find(Vmna(:, kmna) > GlobalConstants.FineTol);
                        if isscalar(vis) && sn.sched(vis) ~= SchedStrategy.INF ...
                                && sn.sched(vis) ~= SchedStrategy.EXT
                            bool = false;
                            reason = sprintf(['The mna method does not support self-looping ' ...
                                'classes (class %d is confined to station %d with no ' ...
                                'inter-station flow to decompose). Use the dec.source method.'], ...
                                kmna, vis);
                            return;
                        end
                    end
                case 'ldqbd'
                    if sn.nclasses ~= 1
                        bool = false;
                        reason = 'The ldqbd method requires a single-class model.';
                        return;
                    end
                case {'inap','inapplus','inapinf','exact'}
                    % RCAT collapses each process to its mean rate: reject
                    % non-exponential; see _kb/06-solver-catalog.md for rationale
                    nonExp = (sn.procid ~= ProcessType.EXP) & isfinite(sn.rates) & (sn.rates > 0);
                    % Exempt signal service entries from the nonExp mask;
                    % see _kb/06-solver-catalog.md for rationale
                    if isfield(sn, 'issignal') && ~isempty(sn.issignal) && any(sn.issignal)
                        issource = false(size(nonExp, 1), 1);
                        for ist = 1:numel(issource)
                            issource(ist) = sn.nodetype(sn.stationToNode(ist)) == NodeType.Source;
                        end
                        nonExp(~issource, logical(sn.issignal(:))') = false;
                    end
                    if any(nonExp(:))
                        [ist, r] = find(nonExp, 1);
                        bool = false;
                        reason = sprintf(['The %s method supports exponential processes only ' ...
                            '(RCAT models each station-class by its mean rate, with no service-phase ' ...
                            'dimension), but station %d class %d is %s. Use the dec.source method ' ...
                            'for non-exponential models.'], method, ist, r, ...
                            ProcessType.toText(sn.procid(ist, r)));
                        return;
                    end
                    % RCAT models every station single-server: reject multi-server;
                    % see _kb/06-solver-catalog.md for rationale
                    multi = isfinite(sn.nservers) & (sn.nservers > 1);
                    if any(multi(:))
                        ist = find(multi, 1);
                        bool = false;
                        reason = sprintf(['The %s method supports single-server stations only ' ...
                            '(RCAT does not model sn.nservers, so a multiserver station is driven ' ...
                            'at rho = lambda/mu instead of lambda/(c*mu)), but station %d has %d ' ...
                            'servers. Use the dec.source method for multiserver models.'], ...
                            method, ist, sn.nservers(ist));
                        return;
                    end
            end
            [bool, reason] = supportsModelMethod@NetworkSolver(self, method);
        end
end
    
    methods (Static)


        function featSupported = getFeatureSet()
            % FEATSUPPORTED = GETFEATURESET()

            featSupported = SolverFeatureSet;
            % MAM features
            featSupported.setTrue({'Sink','Source',...
                'Fork','Join','Forker','Joiner',... % Fork-Join support (via FJ_codes)
                'Delay','DelayStation','Queue',...
                'APH','Coxian','Erlang','Exp','HyperExp','MMPP2','MAP','MMAP','DMAP','ME','RAP',...
                'Det','Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher','InfiniteServer',...
                'ClassSwitch', ...
                'SharedServer','Buffer','Dispatcher',...
                'Server','JobSink','RandomSource','ServiceTunnel',...
                'SchedStrategy_INF','SchedStrategy_PS','SchedStrategy_HOL',...
                'SchedStrategy_FCFSPRPRIO',... % solver_mam_basic: MMAPPH1PRPR
                'SchedStrategy_FCFS',...
                'RoutingStrategy_PROB','RoutingStrategy_RAND',...
                'ClosedClass','SelfLoopingClass',...
                'OpenClass'});
            % Add RCAT (AG) features
            featSupported.setTrue({'Sink', 'Source', ...
                'Fork','Join','Forker','Joiner',... % Fork-Join support (via FJ_codes)
                'Delay', 'DelayStation', 'Queue', ...
                'APH', 'Coxian', 'Erlang', 'Exp', 'HyperExp', ...
                'Det','Gamma','Lognormal','Pareto','Uniform','Weibull',...
                'StatelessClassSwitcher', 'InfiniteServer', ...
                'SharedServer', 'Buffer', 'Dispatcher', ...
                'Server', 'JobSink', 'RandomSource', 'ServiceTunnel', ...
                'SchedStrategy_INF', 'SchedStrategy_PS', ...
                'SchedStrategy_FCFS', ...
                'RoutingStrategy_PROB', 'RoutingStrategy_RAND', ...
                'ClosedClass', ...
                'OpenClass', ...
                'OpenSignal', 'ClosedSignal', ... % G-network signals (solver_mam_ag)
                'SignalType_NEGATIVE', 'SignalType_CATASTROPHE', ...
                'SignalBatchRemoval'}); % AG reads sn.signalremdist
            % Add BMAP/PH/N/N retrial queue features
            featSupported.setTrue({'Retrial', 'BMAP', 'PH'});
            % Setup/delay-off: open stations are solved exactly by
            % qbd_setupdelayoff; closed stations use the per-instance
            % cold-start race of the isfunction branch.
            featSupported.setTrue({'SetupDelayOff'});
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)
            
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverMAM.getFeatureSet();
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
        
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('MAM');
        end

        function libs = getLibrariesUsed(sn, options)
            % GETLIBRARIESUSED Get list of external libraries used by MAM solver
            % Detect libraries used by MAM solver based on topology and method
            libs = {};

            % MAMSolver used for matrix-analytic methods (M/G/1, GI/M/1 types)
            % This includes default and decomposition methods
            if ismember(options.method, {'default', 'dec.source', 'dec.mmap', 'dec.poisson', 'dec.source.mmap'})
                libs{end+1} = 'MAMSolver';
            end

            % Q-MAM used for specific RCAT-based methods
            if ismember(options.method, {'mna', 'inap', 'inapplus', 'inapinf'})
                libs{end+1} = 'Q-MAM';
            end

            % SMCSolver used for QBD (Quasi-Birth-Death) analysis
            % Currently available but not actively used in default paths
            % Uncomment when QBD methods are activated:
            % if ismember(options.method, {'qbd'})
            %     libs{end+1} = 'SMCSolver';
            % end

            % BUTools usage detection (via KPCToolbox MAP functions)
            % BUTools is used when analyzing MAP/PH distributions
            if ~isempty(sn) && isfield(sn, 'proc') && ~isempty(sn.proc) && any(~cellfun(@isempty, sn.proc(:)))
                libs{end+1} = 'BUTools';
            end

            % Remove duplicates and maintain order
            libs = unique(libs, 'stable');
        end
    end
end
