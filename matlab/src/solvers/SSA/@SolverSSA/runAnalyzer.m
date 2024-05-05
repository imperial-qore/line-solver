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

switch options.method
    case 'taussa'
        options.lang = 'java';
end
sn = getStruct(self);

switch options.lang
    case 'java'
        jlineSSASupported = true;
        M = sn.nstations;
        R = sn.nclasses;
        for i=1:M
            for r=1:R
                if ~(sn.procid(i,r) == ProcessType.ID_EXP || sn.procid(i,r) == ProcessType.ID_ERLANG || sn.procid(i,r) == ProcessType.ID_MAP || sn.procid(i,r) == ProcessType.ID_COXIAN || sn.procid(i,r) == ProcessType.ID_DISABLED || sn.procid(i,r) == ProcessType.ID_IMMEDIATE)
                    jlineSSASupported = false;
                end
            end
        end

        if ~jlineSSASupported
            line_error(mfilename,'This model contains features not supported by the Java solver.');
        end

        switch options.method
            case 'default'
                jmodel = LINE2JLINE(self.model);
                jsolver = JLINE.SolverSSA(jmodel);
                import jline.solvers.ssa.*;
                import jline.solvers.ssa.strategies.*;

                M = jmodel.getNumberOfStatefulNodes; %number of stations
                K = jmodel.getNumberOfClasses;    %number of classes

                jsolver.options.samples = options.samples;
                jsolver.options.seed = options.seed;
                tic;
                jsolver.getAvg();
                QN = JLINE.jlinematrix_to_matrix(jsolver.getAvgQLen()); 
                UN = JLINE.jlinematrix_to_matrix(jsolver.getAvgUtil());
                RN = JLINE.jlinematrix_to_matrix(jsolver.getAvgRespT());
                WN = JLINE.jlinematrix_to_matrix(jsolver.getAvgResidT());
                TN = JLINE.jlinematrix_to_matrix(jsolver.getAvgTput());
                AN = JLINE.jlinematrix_to_matrix(jsolver.getAvgArvR());
                CN = JLINE.jlinematrix_to_matrix(jsolver.getAvgSysRespT());
                XN = JLINE.jlinematrix_to_matrix(jsolver.getAvgSysTput());
                runtime = toc;
                self.setAvgResults(QN,UN,RN,TN,AN,WN,CN,XN,runtime);
            case 'taussa'
                line_error(mfilename, 'taussa method is no longer supported.');
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
                line_error(mfilename, 'tauleap method is no longer supported.');
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