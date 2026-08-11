N = 30;
M = 3;
E = 4;
envModel = Environment('MyEnv');
envName = {'Stage1', 'Stage2', 'Stage3', 'Stage4'};
envType = {'UP', 'DOWN', 'FAST' 'SLOW'};

rate = ones(M,E); rate(M,1:E)=(1:E); rate(1,1:E)=(E:-1:1);

qn1 = renv_genqn(rate(:,1),N);
qn2 = renv_genqn(rate(:,2),N);
qn3 = renv_genqn(rate(:,3),N);
qn4 = renv_genqn(rate(:,4),N);
envSubModel = {qn1,qn2,qn3,qn4};
for e=1:E
    envModel.addStage(envName{e}, envType{e}, envSubModel{e});
end

envRates = [0,1,0,0; 0,0,1,1; 1,0,0,1; 1,1,0,0]/2;
for e=1:E
    for h=1:E
        if envRates(e,h)>0
            envModel.addTransition(envName{e}, envName{h}, APH.fitMeanAndSCV(1/envRates(e,h),0.5));
        end
    end
end


%%
fprintf(1,'The metasolver considers an environment with 4 stages and a queueing network with 3 stations.\n')
fprintf(1,'Every time the stage changes, the queueing network will modify the service rates of the stations.\n')

envModel.getStageTable

options = Solver.defaultOptions;
options.timespan = [0,Inf];
options.iter_max = 100;
options.iter_tol = 0.05;
options.method = 'default';

soptions = FLD.defaultOptions;
soptions.timespan = [0,Inf];
soptions.verbose = 0;
envSolver = ENV(envModel, @(model) FLD(model, soptions), options);
[QN,UN,TN] = envSolver.getAvg();
AvgTable = envSolver.getAvgTable()
AvgTables = envSolver.getEnsembleAvgTables()