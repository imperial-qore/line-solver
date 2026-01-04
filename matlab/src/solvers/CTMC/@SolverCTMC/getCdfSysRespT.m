function RD = getCdfSysRespT(self)
% RD = GETCDFRESPT()
%global GlobalConstants.FineTol
sn = self.getStruct;
RD = cell(1, sn.nchains);
N = sn.njobs;
for c=1:sn.nchains
    inchain = sn.inchain{c};
    s = inchain(N(inchain)>0); % tag a class that has non-zero jobs.
    jobclass = self.model.getClassByIndex(s);
    chain = self.model.getClassChain(jobclass);
    [taggedModel, taggedJob] = ModelAdapter.tagChain(chain,jobclass); % diminish jobclass population by 1
    [Q,F,ev] = SolverCTMC(taggedModel,self.options).getGenerator(); % Q: generator, F: filtration, ev: events
    tsn = taggedModel.getStruct;
    tinchain = cell2mat(taggedJob.index);
    D1 = sparse(zeros(size(F{s})));
    for r=tinchain % filter tagged job
        if taggedModel.classes{r}.completes
            for v=1:length(ev)
                if ev{v}.passive{1}.event == EventType.ARV && ev{v}.passive{1}.class == r && ev{v}.passive{1}.node == tsn.refstat(r)
                    D1 = D1 + sparse(F{v});
                end
            end
        end
    end

    D = map_normalize({Q-D1, D1});
    pie_arv = map_pie(D); % state seen upon arrival of a class-r job

    nonZeroRates = abs(Q(Q~=0));
    nonZeroRates = nonZeroRates( nonZeroRates > GlobalConstants.FineTol );
    T = abs(100/min(nonZeroRates)); % solve ode until T = 100 events with the slowest rate
    dT = T/10000; % solve ode until T = 100 events with the slowest rate
    tset = 0:dT:T;

    if year(matlabRelease.Date)>2023 || strcmpi(matlabRelease.Release,"R2023b")
        % use exmpmv
        RD{1,c} = zeros(length(tset),2);
        RD{1,c}(1:length(tset),1) = 1-sum(expmv(D{1}',pie_arv',tset)' ,2);
        RD{1,c}(1:length(tset),2) = tset;
        for t=1:length(tset)
            RD{1,c}(t,1) = 1-pie_arv * expm(D{1}*tset(t)) * ones(length(D{1}),1);
            if RD{1,c}(t,1)>1-GlobalConstants.FineTol
                RD{1,c}(t+1:end,:)=[];
                break
            end
        end
    else
        RD{1,c} = zeros(length(tset),2);
        for t=1:length(tset)
            RD{1,c}(t,2) = tset(t);
            RD{1,c}(t,1) = 1-pie_arv * expm(D{1}*tset(t)) * ones(length(D{1}),1);
            if RD{1,c}(t,1)>1-GlobalConstants.FineTol
                RD{1,c}(t+1:end,:)=[];
                break
            end
        end
    end
end
end