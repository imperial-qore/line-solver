function logData = parseLogs(model,isNodeLogged, metric)
% LOGDATA = PARSELOGS(MODEL,ISNODELOGGED, METRIC)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%line_printf('\nJMT log parsing...');
T0=tic;
sn = model.getStruct;
nclasses = sn.nclasses;
logData = cell(sn.nnodes,sn.nclasses);
nodePreload = sn_get_state_aggr(model.getStruct);
for ind=1:sn.nnodes
    if sn.isstateful(ind) && isNodeLogged(ind)
        logFileArv = [model.getLogPath,sprintf('%s-Arv.csv',model.getNodeNames{ind})];
        logFileDep = [model.getLogPath,sprintf('%s-Dep.csv',model.getNodeNames{ind})];
        %% load arrival process
        if exist(logFileArv,'file') && exist(logFileDep,'file')
            %            logArv=readTable(logFileArv,'Delimiter',';','HeaderLines',1); % raw data
            
            %%
            % unclear if this part works fine if user has another local
            % since JMT might write with a different delimiter
            fid=fopen(logFileArv);
            logArv = textscan(fid, '%s%f%f%s%s%s', 'delimiter',';', 'headerlines',1);
            fclose(fid);
            jobArvTS = logArv{2};
            jobArvID = logArv{3};
            jobArvClass = logArv{4};
            clear logArv;
            
            %%
            jobArvClasses = unique(jobArvClass);
            jobArvClassID = zeros(length(jobArvClass),1);
            for c=1:length(jobArvClasses)
                jobArvClassID(find(strcmp(jobArvClasses{c},jobArvClass))) = findstring(model.getClassNames,jobArvClasses{c});
                %                jobArvClassID(find(strcmp(jobArvClasses{c},jobArvClass)))=c;
            end
            logFileArvMat = [model.getLogPath,filesep,sprintf('%s-Arv.mat',model.getNodeNames{ind})];
            save(logFileArvMat,'jobArvTS','jobArvID','jobArvClass','jobArvClassID');
            
            %% load departure process
            fid=fopen(logFileDep);
            logDep = textscan(fid, '%s%f%f%s%s%s', 'delimiter',';', 'headerlines',1);
            fclose(fid);
            jobDepTS = logDep{2};
            jobDepID = logDep{3};
            jobDepClass = logDep{4};
            clear logDep;
            
            %             jobDepTS = table2array(logDep(:,2));
            %             jobDepID = table2array(logDep(:,3));
            %             jobDepClass = table2cell(logDep(:,4));
            jobDepClasses = unique(jobDepClass);
            jobDepClassID = zeros(length(jobDepClass),1);
            for c=1:length(jobDepClasses)
                jobDepClassID(find(strcmp(jobDepClasses{c},jobDepClass))) = findstring(model.getClassNames,jobDepClasses{c});
                %                jobDepClassID(find(strcmp(jobDepClasses{c},jobDepClass)))=c;
            end
            logFileDepMat = [model.getLogPath,filesep,sprintf('%s-Dep.mat',model.getNodeNames{ind})];
            save(logFileDepMat,'jobDepTS','jobDepID','jobDepClass','jobDepClassID');
            
            nodeState = cell(sn.nnodes,1);
            
            switch metric
                case MetricType.toText(MetricType.QLen)
                    [node_ind, evtype, evclass, evjob] = JMTIO.parseTranState(logFileArvMat, logFileDepMat, nodePreload{ind});
                    
                    %% save in default data structure
                    for r=1:nclasses %0:numOfClasses
                        logData_indr = struct();
                        logData_indr.t = node_ind(:,1);
                        ec = cell(length(evtype),1);
                        for e=1:length(evtype)
                            if evclass(e) == r                                 
                                ec{e,1} = Event(evtype(e), ind, r, NaN, [], node_ind(e,1), evjob(e));
                            else
                                ec{e,1} = [];
                            end
                        end
                        logData_indr.event = ec;
                        logData_indr.QLen = node_ind(:,1+r);
                        %logData_indr.arvID = jobArvID;
                        %logData_indr.depID = jobDepID;                        
                        logData{ind,r} = logData_indr;
                    end
                    nodeState{ind} = node_ind;
                case MetricType.toText(MetricType.RespT)
                    [classResT, jobRespT, jobResTArvTS] = JMTIO.parseTranRespT(logFileArvMat, logFileDepMat);
                    
                    for r=1:nclasses
                        logData{ind,r} = struct();
                        if r <= size(classResT,2)
                            logData{ind,r}.t = jobResTArvTS;
                            logData{ind,r}.RespT = classResT{r};
                            %logData{i,r}.PassT = jobRespT;
                        else
                            logData{ind,r}.t = [];
                            logData{ind,r}.RespT = [];
                            %logData{i,r}.PassT = [];
                        end
                    end
            end
        end
    end
end
runtime=toc(T0);
%line_printf(' completed. Runtime: %f seconds.\n',runtime);
end

