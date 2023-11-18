function model = QN2LINE(sn, modelName)
% MODEL = QN2LINE(QN, MODELNAME)

% Copyright (c) 2012-2023, Imperial College London
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
for i = 1:M
    switch sn.schedid(i)
        case SchedStrategy.ID_INF
            node{i} = DelayStation(model, sn.nodenames{i});
        case SchedStrategy.ID_FORK
            node{i} = ForkStation(model, sn.nodenames{i});
        case SchedStrategy.ID_EXT
            node{i} = Source(model, 'Source'); idSource = i;
            node{M+1} = Sink(model, 'Sink'); hasSink = 1;
        otherwise
            node{i} = Queue(model, sn.nodenames{i}, sn.sched{i});
            node{i}.setNumServers(sn.nservers(i));
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
        for i=1:M
            if sum(nnz(sn.proc{i}{k}{1}))>0
                break
            end
        end
        if isinf(NK(k))
            jobclass{k} = OpenClass(model, sn.classnames{k});
        else
            jobclass{k} = ClosedClass(model, sn.classnames{k}, NK(k), node{i}, 0);
        end
    end
    
    for i=1:M
        SCVik = map_scv(PH{i}{k});
        %        if SCVik >= 0.5
        switch sn.schedid(i)
            case SchedStrategy.ID_EXT
                if isnan(sn.rates(i,k))
                    node{i}.setArrival(jobclass{k}, Disabled.getInstance());
                elseif sn.rates(i,k)==0
                    node{i}.setArrival(jobclass{k}, Immediate.getInstance());
                else
                    node{i}.setArrival(jobclass{k}, APH.fitMeanAndSCV(map_mean(PH{i}{k}),SCVik));
                end
            case SchedStrategy.ID_FORK
                % do nothing
            otherwise
                if isnan(sn.rates(i,k))
                    node{i}.setService(jobclass{k}, Disabled.getInstance());
                elseif sn.rates(i,k)==0
                    node{i}.setService(jobclass{k}, Immediate.getInstance());
                else
                    node{i}.setService(jobclass{k}, APH.fitMeanAndSCV(map_mean(PH{i}{k}),SCVik));
                end
        end
        %        else
        % this could be made more precised by fitting into a 2-state
        % APH, especially if SCV in [0.5,0.1]
        %            nPhases = max(1,round(1/SCVik));
        %            switch sn.schedid(i)
        %                case SchedStrategy.ID_EXT
        %                    node{i}.setArrival(jobclass{k}, Erlang(nPhases/map_mean(PH{i}{k}),nPhases));
        %                case SchedStrategy.ID_FORK
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
        for i=1:M
            for m=1:M
                % routing matrix for each class
                if hasSink && m == idSource % direct to sink
                    myP{k,c}(i,M+1) = rt((i-1)*K+k,(m-1)*K+c);
                else
                    myP{k,c}(i,m) = rt((i-1)*K+k,(m-1)*K+c);
                end
            end
        end
    end
end

model.link(myP);
end
