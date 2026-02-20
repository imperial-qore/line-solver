% This example shows the computation of the generator
N = 2;
M = 2;
E = 3;
envModel = Environment('MyEnv');
envName = {'Stage1', 'Stage2', 'Stage3'};
envType = {'UP', 'DOWN', 'FAST'};

rate = ones(M,E); rate(M,1:E)=(1:E); rate(1,1:E)=(E:-1:1);

qn1 = renv_genqn(rate(:,1),N);
qn2 = renv_genqn(rate(:,2),N);
qn3 = renv_genqn(rate(:,3),N);

envSubModel = {qn1,qn2,qn3};
for e=1:E
    envModel.addStage(envName{e}, envType{e}, envSubModel{e});
end

envRates = circul(3);
for e=1:E
    for h=1:E
        if envRates(e,h)>0
           envModel.addTransition(envName{e}, envName{h}, Erlang.fitMeanAndOrder(1/envRates(e,h),e+h));
        end
    end
end

%%
fprintf(1,'The metasolver considers an environment with 3 stages and a queueing network with 2 stations.\n')
fprintf(1,'This example illustrates the computation of the infinitesimal generator of the system.\n')

envModel.getStageTable

options = Solver.defaultOptions;
options.timespan = [0,Inf];
options.iter_max = 100;
options.iter_tol = 0.05;
options.method = 'default';

soptions = CTMC.defaultOptions;
soptions.timespan = [0,Inf];
soptions.verbose = 0;
soptions.stiff = false;
envSolver = ENV(envModel, @(model) CTMC(model, soptions), options);
[infGen,stageInfGen] = envSolver.getGenerator()

envSolver.getAvgTable()