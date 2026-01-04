clear node jobclass;

figure;
label = {};
nJobs = [1,4,8];
plotColors = {'k';'b--';'r-.'};
for ni = 1:length(nJobs)
    N = nJobs(ni);
    model = Network('model');
    
    node{1} = Delay(model, 'Delay');
    node{2} = Queue(model, 'Queue1', SchedStrategy.PS);
    node{3} = Queue(model, 'Queue2', SchedStrategy.PS);
    
    jobclass{1} = ClosedClass(model, 'Class1', N, node{1}, 0);
    
    jobclass{1}.completes = false;
    
    node{1}.setService(jobclass{1}, Exp(1/1));
    node{2}.setService(jobclass{1}, Exp(1/2));
    node{3}.setService(jobclass{1}, Exp(1/2));
    
    M = model.getNumberOfStations();
    K = model.getNumberOfClasses();
    
    P = circul(M);
    model.link(P);
    %%
    options = FLD.defaultOptions;
    options.iter_max = 100;
    solver = FLD(model, options);
    AvgRespT{ni} = solver.getAvgRespT
    RDfluid = solver.getCdfRespT();
    %%
    for c=1:model.getNumberOfClasses
        for i=1:model.getNumberOfStations
            AvgRespTfromCDF{ni}(i,c) = diff(RDfluid{i,c}(:,1))'*RDfluid{i,c}(2:end,2); %mean
            PowerMoment2_R{ni}(i,c) = diff(RDfluid{i,c}(:,1))'*(RDfluid{i,c}(2:end,2).^2);
            Variance_R{ni}(i,c) = PowerMoment2_R{ni}(i,c)-AvgRespTfromCDF{ni}(i,c)^2; %variance
            SqCoeffOfVariationRespTfromCDF{ni}(i,c) = (Variance_R{ni}(i,c))/AvgRespTfromCDF{ni}(i,c)^2; %scv
        end
    end
    for i=2
        for c=1:model.getNumberOfClasses
            semilogx(RDfluid{i,c}(:,2),RDfluid{i,c}(:,1),plotColors{find(N==nJobs)}); hold on;
        end
    end
    AvgRespTfromCDF{ni}
    SqCoeffOfVariationRespTfromCDF{ni}
    label{end+1} = ['N=', num2str(N),' jobs'];    
end
legend(label,'Location','SouthEast');
xlim([0,200])
title('Response time CDF at station 3 under increasing populations')
ylabel('Pr[RespT < t]')
xlabel('Response time t')
