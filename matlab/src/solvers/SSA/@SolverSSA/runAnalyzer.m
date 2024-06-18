function [runtime, tranSysState, tranSync] = runAnalyzer(self, options)
% [RUNTIME, TRANSYSSTATE] = RUN()

T0=tic;
if nargin<2 %~exist('options','var')
    options = self.getOptions;
end

if self.enableChecks && ~self.supports(self.model)
    line_error(mfilename,'This model contains features not supported by the solver.');
end

self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);

%options.lang = 'java';

switch options.method
    case 'taussa'
        options.lang = 'java';
    case 'default'
        switch options.lang
            case 'java'
                options.verbose = VerboseLevel.SILENT;
                options.method = 'para';
            case 'matlab'
                options.method = 'serial';
        end
end
sn = getStruct(self);

switch options.lang
    case 'java'
        switch options.method
            case {'default','serial','para','parallel'}
                jmodel = LINE2JLINE(self.model);
                joptions = jline.solvers.SolverOptions(jline.lang.constant.SolverType.SSA);
                joptions.samples = options.samples;
                joptions.seed = options.seed;
                tic;
                jsolver = JLINE.SolverSSA(jmodel, options);
                [QN,UN,RN,WN,AN,TN] = JLINE.arrayListToResults(jsolver.getAvgTable);           
                CN = JLINE.jlinematrix_to_matrix(jsolver.getAvgSysRespT());
                XN = JLINE.jlinematrix_to_matrix(jsolver.getAvgSysTput());
                runtime = toc;
                self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime);
            case 'taussa'
                line_error(mfilename, 'taussa method is temporarily disabled.');
                % [XN,UN,QN,RN,TN,CN, tranSysState, tranSync] = solver_ssa_analyzer_taussa(self.model, options, 0, 0.0);
                % runtime = toc(T0);
                % T = getAvgTputHandles(self);
                % if ~isempty(T) && ~isempty(TN)
                %     AN = zeros(M,R);
                %     for i=1:M
                %         for j=1:M
                %             for k=1:R
                %                 for r=1:R
                %                     AN(i,k) = AN(i,k) + TN(j,r)*sn.rt((j-1)*R+r, (i-1)*R+k);
                %                 end
                %             end
                %         end
                %     end
                % else
                %     AN = [];
                % end
                % self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime);
            case {'tauleap'}
                line_error(mfilename, 'tauleap method is temporarily disabled.');
                % [XN,UN,QN,RN,TN,CN, tranSysState, tranSync] = solver_ssa_analyzer_taussa(self.model, options, 1, 0.0);
                % runtime = toc(T0);
                % sn = self.getStruct;
                % M = sn.nstations;
                % R = sn.nclasses;
                % T = getAvgTputHandles(self);
                % if ~isempty(T) && ~isempty(TN)
                %     AN = zeros(M,R);
                %     for i=1:M
                %         for j=1:M
                %             for k=1:R
                %                 for r=1:R
                %                     AN(i,k) = AN(i,k) + TN(j,r)*sn.rt((j-1)*R+r, (i-1)*R+k);
                %                 end
                %             end
                %         end
                %     end
                % else
                %     AN = [];
                % end
                % self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime);
        end
    case 'matlab'
        [QN,UN,RN,TN,CN,XN,~, tranSysState, tranSync, sn] = solver_ssa_analyzer(sn, options);

        for isf=1:sn.nstateful
            ind = sn.statefulToNode(isf);
            switch sn.nodetype(sn.statefulToNode(isf))
                case NodeType.Cache
                    self.model.nodes{sn.statefulToNode(isf)}.setResultHitProb(sn.nodeparam{ind}.actualhitprob);
                    self.model.nodes{sn.statefulToNode(isf)}.setResultMissProb(sn.nodeparam{ind}.actualmissprob);
                    self.model.refreshChains();
            end
        end
        runtime = toc(T0);
        M = sn.nstations;
        R = sn.nclasses;
        T = getAvgTputHandles(self);
        if ~isempty(T) && ~isempty(TN)
            AN = zeros(M,R);
            for i=1:M
                for j=1:M
                    for k=1:R
                        for r=1:R
                            AN(i,k) = AN(i,k) + TN(j,r)*sn.rt((j-1)*R+r, (i-1)*R+k);
                        end
                    end
                end
            end
        else
            AN = [];
        end
        self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime);
        self.result.space = sn.space;
end
end