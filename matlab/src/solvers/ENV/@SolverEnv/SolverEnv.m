classdef SolverEnv < EnsembleSolver
    % Solver for models immersed in a random environment.
    %
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.


    properties
        env; % user-supplied representation of each stage transition
        envObj;
        sn;
        resetFromMarginal;
    end

    methods
        function self = SolverEnv(renv, solverFactory, options)
            % SELF = SOLVERENV(ENV,SOLVERFACTORY,OPTIONS)
            self@EnsembleSolver(renv, mfilename);
            if nargin>=3 %exist('options','var')
                self.setOptions(options);
            else
                self.setOptions(SolverEnv.defaultOptions);
            end            
            self.envObj = renv;
            self.ensemble = renv.getEnsemble;
            env = renv.getEnv;
            self.env = env;
            for e=1:length(self.env)
                self.sn{e} = self.ensemble{e}.getStruct;
                self.setSolver(solverFactory(self.ensemble{e}),e);
            end

            for e=1:length(self.env)
                for h=1:length(self.env)
                    self.resetFromMarginal{e,h} = renv.resetFun{e,h};
                end
            end

            for e=1:length(self.env)
                for h=1:length(self.env)
                    if isa(self.env{e,h},'Disabled')
                        self.env{e,h} = Exp(0);
                    elseif ~isa(self.env{e,h},'MarkovianDistribution')
                        line_error(mfilename,sprintf('The distribution of the environment transition from stage %d to %d is not supported by the %s solver.',e,h,self.getName));
                    end
                end
            end

            for e=1:length(self.ensemble)
                if ~self.solvers{e}.supports(self.ensemble{e})
                    line_error(mfilename,sprintf('Model in the environment stage %d is not supported by the %s solver.',e,self.getName));
                end
            end
        end

        function bool = converged(self, it) % convergence test at iteration it
            % BOOL = CONVERGED(IT) % CONVERGENCE TEST AT ITERATION IT

            bool = false;
            E = self.getNumberOfModels;
            if it>1
                mapes = zeros(1,E);
                for e=1:E
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for j=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            % error is calculated only on entry value (t=0)
                            if ~isempty(self.results{it,e}.Tran.Avg.Q{i,j})
                                mapes(e) = max(mapes(e), mape(self.results{it,e}.Tran.Avg.Q{i,j}.metric(1), self.results{it-1,e}.Tran.Avg.Q{i,j}.metric(1)));
                            end
                        end
                    end
                end
                if max(mapes) < self.options.iter_tol
                    bool = true;
                end
            end
        end


        function init(self)
            % INIT()
            options = self.options;
            if isfield(options,'seed')
                Solver.resetRandomGeneratorSeed(options.seed);
            end
            self.envObj.init();
        end

        function pre(self, it)
            % PRE(IT)

            if it==1
                for e=list(self)
                    if isinf(self.getSolver(e).options.timespan(2))
                        [QN,~,~,~] = self.getSolver(e).getAvg();
                    else
                        [QNt,~,~] = self.getSolver(e).getTranAvg();
                        QN = cellfun(@(q) q.metric(end), QNt);
                    end
                    self.ensemble{e}.initFromMarginal(QN);
                end
            end
        end

        % solves model in stage e
        function [results_e, runtime] = analyze(self, it, e)
            % [RESULTS_E, RUNTIME] = ANALYZE(IT, E)

            results_e = struct();
            results_e.('Tran') = struct();
            results_e.Tran.('Avg') = [];
            T0 = tic;
            runtime = toc(T0);
            %% initialize
            [Qt,Ut,Tt] = self.ensemble{e}.getTranHandles;
            %[results_e.Avg.Q, results_e.Avg.U, results_e.Avg.R, results_e.Avg.T] = self.solvers{e}.getAvg();
            self.solvers{e}.resetResults;
            [QNt,UNt,TNt] = self.solvers{e}.getTranAvg(Qt,Ut,Tt);
            results_e.Tran.Avg.Q = QNt; %cellfun(@(c) c.metric, QNt,'UniformOutput',false);
            results_e.Tran.Avg.U = UNt; %cellfun(@(c) c.metric, UNt,'UniformOutput',false);
            results_e.Tran.Avg.T = TNt; %cellfun(@(c) c.metric, TNt,'UniformOutput',false);
            %[results_e.Tran.Avg.Q, results_e.Tran.Avg.U, results_e.Tran.Avg.T] = self.solvers{e}.getTranAvg(Qt,Ut,Tt);
        end

        function post(self, it)
            % POST(IT)

            E = self.getNumberOfModels;
            for e=1:E
                for h = 1:E
                    Qexit{e,h} = zeros(size(self.results{it,e}.Tran.Avg.Q));
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            if ~isempty(self.results{it,e}.Tran.Avg.Q{i,r})
                                w{e,h} = [0, map_cdf(self.envObj.proc{e}{h}, self.results{it,e}.Tran.Avg.Q{i,r}.t(2:end)) - map_cdf(self.envObj.proc{e}{h}, self.results{it,e}.Tran.Avg.Q{i,r}.t(1:end-1))]';
                                if ~isnan(w{e,h})
                                    Qexit{e,h}(i,r) = self.results{it,e}.Tran.Avg.Q{i,r}.metric'*w{e,h}/sum(w{e,h});
                                else
                                    Qexit{e,h}(i,r) = 0;
                                end
                            else
                                w{e,h} = 0;
                            end
                        end
                    end
                end
            end

            Qentry = cell(1,E); % average entry queue-length
            for e = 1:E
                Qentry{e} = zeros(size(Qexit{e}));
                for h=1:E
                    % probability of coming from h to e \times resetFun(Qexit from h to e
                    if self.envObj.probOrig(h,e) > 0
                        Qentry{e} = Qentry{e} + self.envObj.probOrig(h,e) * self.resetFromMarginal{h,e}(Qexit{h,e});
                    end
                end
                self.solvers{e}.resetResults();
                self.ensemble{e}.initFromMarginal(Qentry{e});
            end
        end

        function finish(self)
            % FINISH()

            it = size(self.results,1); % use last iteration
            E = self.getNumberOfModels;
            for e=1:E
                QExit{e}=[];
                UExit{e}=[];
                TExit{e}=[];
                if it>0
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            if ~isempty(self.results{it,e}.Tran.Avg.Q{i,r})
                                w{e} = [0, map_cdf(self.envObj.holdTime{e}, self.results{it,e}.Tran.Avg.Q{i,r}.t(2:end)) - map_cdf(self.envObj.holdTime{e}, self.results{it,e}.Tran.Avg.Q{i,r}.t(1:end-1))]';
                                QExit{e}(i,r) = self.results{it,e}.Tran.Avg.Q{i,r}.metric'*w{e}/sum(w{e});
                                UExit{e}(i,r) = self.results{it,e}.Tran.Avg.U{i,r}.metric'*w{e}/sum(w{e});
                                TExit{e}(i,r) = self.results{it,e}.Tran.Avg.T{i,r}.metric'*w{e}/sum(w{e});
                            else
                                QExit{e}(i,r) = 0;
                                UExit{e}(i,r) = 0;
                                TExit{e}(i,r) = 0;
                            end
                        end
                    end
                end
                %                 for h = 1:E
                %                     QE{e,h} = zeros(size(self.results{it,e}.Tran.Avg.Q));
                %                     for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                %                         for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                %                             w{e,h} = [0, map_cdf(self.envObj.proc{e}{h}, self.results{it,e}.Tran.Avg.Q{i,r}(2:end,2)) - map_cdf(self.envObj.proc{e}{h}, self.results{it,e}.Tran.Avg.Q{i,r}(1:end-1,2))]';
                %                             if ~isnan(w{e,h})
                %                                 QE{e,h}(i,r) = self.results{it,e}.Tran.Avg.Q{i,r}(:,1)'*w{e,h}/sum(w{e,h});
                %                             else
                %                                 QE{e,h}(i,r) = 0;
                %                             end
                %                         end
                %                     end
                %                 end
            end

            Qval=0*QExit{e};
            Uval=0*UExit{e};
            Tval=0*TExit{e};
            for e=1:E
                Qval = Qval + self.envObj.probEnv(e) * QExit{e}; % to check
                Uval = Uval + self.envObj.probEnv(e) * UExit{e}; % to check
                Tval = Tval + self.envObj.probEnv(e) * TExit{e}; % to check
            end
            self.result.Avg.Q = Qval;
            %    self.result.Avg.R = R;
            %    self.result.Avg.X = X;
            self.result.Avg.U = Uval;
            self.result.Avg.T = Tval;
            %    self.result.Avg.C = C;
            %self.result.runtime = runtime;
            %if self.options.verbose
            %    line_printf('\n');
            %end
        end

        function name = getName(self)
            % NAME = GETNAME()

            name = mfilename;
        end

        function [renvInfGen, stageInfGen, renvEventFilt, stageEventFilt, renvEvents, stageEvents] = getGenerator(self)
            % [renvInfGen, stageInfGen, stageEventFilt, stageEvents] = getGenerator(self)

            E = self.getNumberOfModels;
            stageInfGen = cell(1,E);
            stageEventFilt = cell(1,E);
            stageEvents = cell(1,E);
            for e=1:E
                if isa(self.solvers{e},'SolverCTMC')
                    [stageInfGen{e}, stageEventFilt{e}, stageEvents{e}] = self.solvers{e}.getGenerator();
                else
                    line_error(mfilename,'This method requires SolverENV to be instantiated with the CTMC solver.');
                end
            end

            nstates = cellfun(@length, stageInfGen);
            nphases = cellfun(@(p) p.getNumberOfPhases, self.env);
            nphases = nphases - eye(length(nphases));

            renvInfGen = cell(E,E);

            for e=1:E
                renvInfGen{e,e} = (stageInfGen{e});
                for h=1:E
                    if h~=e
                        resetMatrix_eh = eye(nstates(e), nstates(h));
                        renvInfGen{e,h} = resetMatrix_eh;
                    end
                end
            end

            renvEvents = cell(1,0);
            for e=1:E
                for h=1:E
                    if h~=e
                        renvInfGen{e,e} = krons(renvInfGen{e,e}, (self.env{e,h}.getPH{1}));
                        pie = map_pie(self.env{h,e}.getPH);
                        if isnan(pie)
                            pie=ones(size(pie)); % kron so 1 means ignore
                        end
                        renvInfGen{e,h} = kron(renvInfGen{e,h}, sparse((self.env{e,h}.getPH{2}) * ones(nphases(e,h),1) * pie));
                        for i=1:self.ensemble{e}.getNumberOfNodes
                            renvEvents{1,end+1} = Event(EventType.ID_STAGE, i, NaN, NaN, [e,h]);  %#ok<AGROW>
                        end
                        for f=1:E
                            if f~=h && f~=e
                                pie = map_pie(self.env{f,h}.getPH);
                                if isnan(pie)
                                    pie=ones(size(pie)); % kron so 1 means ignore
                                end
                                renvInfGen{e,f} = kron(renvInfGen{e,f}, ones(nphases(e,h),1) * pie);
                            end
                        end
                    end
                end
            end

            if nargout>2
                line_warning(mfilename, 'Some of the requested output parameters are not implemented yet.\n');
            end

            renvEventFilt = cell(E,E);
            for e=1:E
                for h=1:E
                    tmpCell = renvInfGen;
                    for e1=1:E
                        tmpCell{e1,e1} = tmpCell{e1,e1} * 0;
                        for h1=1:E
                            if e~= e1 && h~=h1
                                tmpCell{e1,h1} = tmpCell{e1,h1} * 0;
                            end
                        end
                    end
                    renvEventFilt{e,h} = cell2mat(tmpCell);
                end
            end

            renvInfGen = cell2mat(renvInfGen);
            renvInfGen = ctmc_makeinfgen(renvInfGen);
        end

        function varargout = getAvg(varargin)
            % [QNCLASS, UNCLASS, TNCLASS] = GETAVG()
            [varargout{1:nargout}] = getEnsembleAvg( varargin{:} );
        end

        function [QNclass, UNclass, RNclass, TNclass, ANclass, WNclass] = getEnsembleAvg(self)
            % [QNCLASS, UNCLASS, TNCLASS] = GETENSEMBLEAVG()
            RNclass=[];
            ANclass=[];
            WNclass=[];

            if isempty(self.result) || (isfield(self.options,'force') && self.options.force)
                iterate(self);
                if isempty(self.result)
                    QNclass=[];
                    UNclass=[];
                    TNclass=[];
                    return
                end
            end
            QNclass = self.result.Avg.Q;
            UNclass = self.result.Avg.U;
            TNclass = self.result.Avg.T;
            WNclass = QNclass ./ TNclass;
            RNclass = NaN*WNclass;
            ANclass = NaN*TNclass;
        end

        function [AvgTable,QT,UT,TT] = getAvgTable(self,keepDisabled)
            % [AVGTABLE,QT,UT,TT] = GETAVGTABLE(SELF,KEEPDISABLED)
            % Return table of average station metrics

            if nargin<2 %if ~exist('keepDisabled','var')
                keepDisabled = false;
            end

            [QN,UN,~,TN] = getAvg(self);
            M = size(QN,1);
            K = size(QN,2);
            Q = self.result.Avg.Q;
            U = self.result.Avg.U;
            T = self.result.Avg.T;
            if isempty(QN)
                AvgTable = Table();
                QT = Table();
                UT = Table();
                TT = Table();
            elseif ~keepDisabled
                Qval = []; Uval = []; Tval = [];
                JobClass = {};
                Station = {};
                for i=1:M
                    for k=1:K
                        if any(sum([QN(i,k),UN(i,k),TN(i,k)])>0)
                            JobClass{end+1,1} = self.sn{1}.classnames{k};
                            Station{end+1,1} = self.sn{1}.nodenames{self.sn{1}.stationToNode(i)};
                            Qval(end+1) = QN(i,k);
                            Uval(end+1) = UN(i,k);
                            Tval(end+1) = TN(i,k);
                        end
                    end
                end
                QLen = Qval(:); % we need to save first in a variable named like the column
                QT = Table(Station,JobClass,QLen);
                Util = Uval(:); % we need to save first in a variable named like the column
                UT = Table(Station,JobClass,Util);
                Tput = Tval(:); % we need to save first in a variable named like the column
                Station = categorical(Station);
                JobClass = categorical(JobClass);
                TT = Table(Station,JobClass,Tput);
                RespT = QLen ./ Tput;
                AvgTable = Table(Station,JobClass,QLen,Util,RespT,Tput);
            else
                Qval = zeros(M,K); Uval = zeros(M,K);
                JobClass = cell(K*M,1);
                Station = cell(K*M,1);
                for i=1:M
                    for k=1:K
                        JobClass{(i-1)*K+k} = Q{i,k}.class.name;
                        Station{(i-1)*K+k} = Q{i,k}.station.name;
                        Qval((i-1)*K+k) = QN(i,k);
                        Uval((i-1)*K+k) = UN(i,k);
                        Tval((i-1)*K+k) = TN(i,k);
                    end
                end
                Station = categorical(Station);
                JobClass = categorical(JobClass);
                QLen = Qval(:); % we need to save first in a variable named like the column
                QT = Table(Station,JobClass,QLen);
                Util = Uval(:); % we need to save first in a variable named like the column
                UT = Table(Station,JobClass,Util);
                Tput = Tval(:); % we need to save first in a variable named like the column
                TT = Table(Station,JobClass,Tput);
                RespT = QLen ./ Tput;
                AvgTable = Table(Station,JobClass,QLen,Util,RespT,Tput);
            end
        end
    end

    methods (Static)
        function [allMethods] = listValidMethods()
            % allMethods = LISTVALIDMETHODS()
            % List valid methods for this solver
            allMethods = {'default'};
        end
        
        function [bool, featSupported] = supports(model)
            % [BOOL, FEATSUPPORTED] = SUPPORTS(MODEL)

            featUsed = model.getUsedLangFeatures();

            featSupported = SolverFeatureSet;

            % Nodes
            featSupported.setTrue('ClassSwitch');
            featSupported.setTrue('DelayStation');
            featSupported.setTrue('Queue');
            featSupported.setTrue('Sink');
            featSupported.setTrue('Source');

            % Distributions
            featSupported.setTrue('Coxian');
            featSupported.setTrue('Cox2');
            featSupported.setTrue('Erlang');
            featSupported.setTrue('Exp');
            featSupported.setTrue('HyperExp');

            % Sections
            featSupported.setTrue('StatelessClassSwitcher'); % Section
            featSupported.setTrue('InfiniteServer'); % Section
            featSupported.setTrue('SharedServer'); % Section
            featSupported.setTrue('Buffer'); % Section
            featSupported.setTrue('Dispatcher'); % Section
            featSupported.setTrue('Server'); % Section (Non-preemptive)
            featSupported.setTrue('JobSink'); % Section
            featSupported.setTrue('RandomSource'); % Section
            featSupported.setTrue('ServiceTunnel'); % Section

            % Scheduling strategy
            featSupported.setTrue('SchedStrategy_INF');
            featSupported.setTrue('SchedStrategy_PS');
            featSupported.setTrue('RoutingStrategy_PROB');
            featSupported.setTrue('RoutingStrategy_RAND');
            featSupported.setTrue('RoutingStrategy_RROBIN'); % with SolverJMT

            % Customer Classes
            featSupported.setTrue('ClosedClass');
            featSupported.setTrue('OpenClass');

            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
    end

    methods (Static)
        % ensemble solver options
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            options = lineDefaults('Env');
        end
    end    
end

