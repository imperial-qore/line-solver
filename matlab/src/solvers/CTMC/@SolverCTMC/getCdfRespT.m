function RD = getCdfRespT(self, R)
% RD = GETCDFRESPT(R)

if GlobalConstants.DummyMode
    RD = cell(1,R);
    return
end

if nargin<2 %~exist('R','var')
    R = getAvgRespTHandles(self);
end
sn = self.getStruct;
RD = cell(sn.nstations, sn.nclasses);
M = sn.nstations;
K = sn.nclasses;
N = sn.njobs;
if any(isinf(N))
    line_error(mfilename,'getCdfRespT is presently supported only for closed models.')
end
for c=1:sn.nchains
    inchain = sn.inchain{c};
    s = inchain(N(inchain)>0); % tag a class that has non-zero jobs.
    jobclass = self.model.getClassByIndex(s);
    chain = self.model.getClassChain(jobclass);
    [taggedModel, taggedJob] = ModelAdapter.tagChain(chain,jobclass); % diminish jobclass population by 1
    %taggedModel.stations{:}
    [Q,F,ev] = SolverCTMC(taggedModel,self.options).getGenerator(); % Q: generator, F: filtration, ev: events
    tsn = taggedModel.getStruct;
    tinchain = cell2mat(taggedJob.index);

    for ist=1:M
        for ir=1:length(tinchain)
            r = tinchain(ir);
            A1 = sparse(zeros(size(F{s})));
            for r=tinchain % filter tagged job
                if taggedModel.classes{r}.completes
                    for v=1:length(ev)
                        if ev{v}.passive{1}.event == EventType.ARV && ev{v}.passive{1}.class == r && ev{v}.passive{1}.node == ist
                            A1 = A1 + sparse(F{v});
                        end
                    end
                end
            end

            D1 = sparse(zeros(size(F{s})));
            for r=tinchain % filter tagged job
                if taggedModel.classes{r}.completes
                    for v=1:length(ev)
                        if ev{v}.active{1}.event == EventType.DEP && ev{v}.active{1}.class == r && ev{v}.active{1}.node == ist
                            D1 = D1 + sparse(F{v});
                        end
                    end
                end
            end

            A = map_normalize({Q-A1, A1});
            pie_arv = map_pie(A); % state seen upon arrival of a class-r job
            D = map_normalize({Q-D1, D1});

            nonZeroRates = abs(Q(Q~=0));
            nonZeroRates = nonZeroRates( nonZeroRates > GlobalConstants.FineTol );
            T = abs(100/min(nonZeroRates)); % solve ode until T = 100 events with the slowest rate
            dT = T/100000; % solve ode until T = 100 events with the slowest rate
            tset = 0:dT:T;

            rorig = chain.index{ir};
            RD{ist,rorig} = zeros(length(tset),2);
            tic
            if year(matlabRelease.Date)>2023 || strcmpi(matlabRelease.Release,"R2023b")
                % use exmpmv
                RD{ist,rorig}(1:length(tset),1) = 1-sum(expmv(D{1}',pie_arv',tset)' ,2);
                RD{ist,rorig}(1:length(tset),2) = tset;
                for t=1:length(tset)
                    if RD{ist,rorig}(t,1)>1-GlobalConstants.CoarseTol
                        RD{ist,rorig}(t+1:end,:)=[];
                        break
                    end
                end
            else % no expmv function available
                for t=1:length(tset)
                    RD{ist,rorig}(t,2) = tset(t);
                    RD{ist,rorig}(t,1) = 1-pie_arv * expm(D{1}*tset(t)) * ones(length(D{1}),1);
                    if RD{ist,rorig}(t,1)>1-GlobalConstants.CoarseTol
                        RD{ist,rorig}(t+1:end,:)=[];
                        break
                    end
                end
            end
            toc
        end
    end
end
end