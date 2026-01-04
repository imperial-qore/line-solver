function model = QN2LINE(sn, modelName)
% MODEL = QN2LINE(QN, MODELNAME)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin<2%~exist('modelName','var')
    modelName = 'sn';
end
%%
M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
rt = sn.rt;
NK = sn.njobs;  % initial population per class
Ktrue = nnz(NK); % classes that are not artificial

%% initialization
model = Network(modelName);
hasSink = 0;
idSource = [];
for ist = 1:M
    switch sn.sched(ist)
        case SchedStrategy.INF
            node{ist} = Delay(model, sn.nodenames{ist});
        case SchedStrategy.FORK
            node{ist} = ForkStation(model, sn.nodenames{ist});
        case SchedStrategy.EXT
            node{ist} = Source(model, 'Source'); idSource = ist;
            node{M+1} = Sink(model, 'Sink'); hasSink = 1;
        otherwise
            node{ist} = Queue(model, sn.nodenames{ist}, SchedStrategy.toFeature(sn.sched(ist)));
            node{ist}.setNumServers(sn.nservers(ist));
    end
end

PH = sn.proc;
for k = 1:K
    if k<=Ktrue
        if isinf(NK(k))
            jobclass{k} = OpenClass(model, sn.classnames{k}, 0);
        else
            jobclass{k} = ClosedClass(model, sn.classnames{k}, NK(k), node{sn.refstat(k)}, 0);
        end
    else
        % if the reference node is unspecified, as in artificial classes,
        % set it to the first node where the rate for this class is
        % non-null
        for ist=1:M
            if sum(nnz(sn.proc{ist}{k}{1}))>0
                break
            end
        end
        if isinf(NK(k))
            jobclass{k} = OpenClass(model, sn.classnames{k});
        else
            jobclass{k} = ClosedClass(model, sn.classnames{k}, NK(k), node{ist}, 0);
        end
    end
    
    for ist=1:M
        SCVik = map_scv(PH{ist}{k});
        %        if SCVik >= 0.5
        switch sn.sched(ist)
            case SchedStrategy.EXT
                if isnan(sn.rates(ist,k))
                    node{ist}.setArrival(jobclass{k}, Disabled.getInstance());
                elseif sn.rates(ist,k)==0
                    node{ist}.setArrival(jobclass{k}, Immediate.getInstance());
                else
                    node{ist}.setArrival(jobclass{k}, APH.fitMeanAndSCV(map_mean(PH{ist}{k}),SCVik));
                end
            case SchedStrategy.FORK
                % do nothing
            otherwise
                if isnan(sn.rates(ist,k))
                    node{ist}.setService(jobclass{k}, Disabled.getInstance());
                elseif sn.rates(ist,k)==0
                    node{ist}.setService(jobclass{k}, Immediate.getInstance());
                else
                    node{ist}.setService(jobclass{k}, APH.fitMeanAndSCV(map_mean(PH{ist}{k}),SCVik));
                end
        end
        %        else
        % this could be made more precised by fitting into a 2-state
        % APH, especially if SCV in [0.5,0.1]
        %            nPhases = max(1,round(1/SCVik));
        %            switch sn.sched(i)
        %                case SchedStrategy.EXT
        %                    node{i}.setArrival(jobclass{k}, Erlang(nPhases/map_mean(PH{i}{k}),nPhases));
        %                case SchedStrategy.FORK
        % do nothing
        %                otherwise
        %                    node{i}.setService(jobclass{k}, Erlang(nPhases/map_mean(PH{i}{k}),nPhases));
        %            end
    end
    %end
end

myP = cell(K,K);
for k = 1:K
    for c = 1:K
        myP{k,c} = zeros(M+hasSink);
        for ist=1:M
            for m=1:M
                % routing matrix for each class
                if hasSink && m == idSource % direct to sink
                    myP{k,c}(ist,M+1) = rt((ist-1)*K+k,(m-1)*K+c);
                else
                    myP{k,c}(ist,m) = rt((ist-1)*K+k,(m-1)*K+c);
                end
            end
        end
    end
end

model.link(myP);
end
