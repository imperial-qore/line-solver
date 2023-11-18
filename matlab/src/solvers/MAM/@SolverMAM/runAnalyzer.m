function runtime = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver

T0=tic;
if nargin<2
    options = self.getOptions;
end

if self.enableChecks && ~self.supports(self.model)
    line_warning(mfilename,'This model contains features not supported by the solver.\n');
    line_error(mfilename,'This model contains features not supported by the solver.\n');
end

self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);

switch options.lang
    case 'java'
        jmodel = LINE2JLINE(self.model);
        M = jmodel.getNumberOfStatefulNodes;
        R = jmodel.getNumberOfClasses;
        jsolver = JLINE.SolverMAM(jmodel);
        [QN,UN,RN,~,TN] = JLINE.arrayListToResults(jsolver.getAvgTable);
        runtime = toc(T0);
        CN = [];
        XN = [];
        QN = reshape(QN',R,M)';
        UN = reshape(UN',R,M)';
        RN = reshape(RN',R,M)';
        TN = reshape(TN',R,M)';
        lG = NaN;
        lastiter = NaN;
        sn = self.getStruct;
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
        self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,options.method,lastiter);
        self.result.Prob.logNormConstAggr = lG;
        return
    case 'matlab'
        sn = getStruct(self);

        [QN,UN,RN,TN,CN,XN] = solver_mam_analyzer(sn, options);
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
        runtime=toc(T0);
        self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime);
end
end