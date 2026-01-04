function [classResT, jobRespT, jobResTArvTS] = parseTranRespT(fileArv, fileDep)
% [CLASSREST, JOBRESPT, JOBRESTARVTS] = PARSETRANRESPT(FILEARV, FILEDEP)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

load(fileArv,'jobArvTS','jobArvClassID','jobArvID');
load(fileDep,'jobDepTS','jobDepClassID','jobDepID');
%% compute job residence times at the station
jobIDs = unique([jobArvID;jobDepID]);
% We count the total number of observed visits. Jobs not departed do not
% complete the visit and are therefore not counted.
jobRespT = cell(length(jobIDs),1);
jobArvClass = cell(length(jobIDs),1);
jobResTArvTS = cell(length(jobIDs),1); % arrival timestamp associated to jobResT entry
jobLogData = [jobDepTS, jobDepID, -ones(size(jobDepTS)), jobDepClassID; jobArvTS, jobArvID, ones(size(jobArvTS)), jobArvClassID];
jobLogData = sortrows(jobLogData,2); % sorted by job ID

ptrJob=zeros(length(jobIDs),1);
ptrJob(1) = 1;
for j=2:size(jobLogData,1)
    if jobLogData(j,2)>jobLogData(j-1,2) % change of job ID
        ptrJob(1+jobLogData(j,2))=j;
    end
end
ptrJob(end+1)=length(jobLogData)+1;

uIDs=unique(jobIDs)';
for id=uIDs
    jobData = jobLogData(ptrJob(1+id):(ptrJob(1+id+1)-1),[1,3]);
    jobClassData{1+id} = jobLogData(ptrJob(1+id):(ptrJob(1+id+1)-1),4);
    %if ~isempty(jobClassData{1+id})
    try
        if length(jobData)>2
            jobData = sortrows(jobData,1);
            firstTS=find(jobData(:,2)>0,1,'first'); % ts of first arrival
            lastTS=find(jobData(:,2)<0,1,'last'); % ts of last departure
            jobData=jobData(firstTS:lastTS,1);
            jobClassData{1+id}=jobClassData{1+id}(firstTS:lastTS,1);
            if length(jobData(2:2:end)) == length(jobData(1:2:end))
                jobResTArvTS{1+id}=jobData(1:2:end); % unclear what happens in case of dropping
                jobRespT{1+id}=jobData(2:2:end)-jobData(1:2:end); % unclear what happens in case of dropping
                jobArvClass{1+id} = jobClassData{1+id}(1:2:end);
            else
                jobResTArvTS{1+id}=jobData(1:2:end-1); % unclear what happens in case of dropping
                jobRespT{1+id}=jobData(2:2:end-1)-jobData(1:2:end-1); % unclear what happens in case of dropping
                jobArvClass{1+id} = jobClassData{1+id}(1:2:end-1);
            end
        else
            jobResTArvTS{1+id}=min(jobData(:,1));
            jobRespT{1+id}=max(jobData(:,1))-min(jobData(:,1));
            jobArvClass{1+id} = jobClassData{1+id}(1);
        end
    %catch
    %    keyboard
    end
    %end
end

%% compute per-class residence times at the station
ID2Class = unique([jobDepID,jobDepClassID],'rows');
for c=unique(ID2Class(:,2))'
    classResT{c} = [];
    jobsInClass = ID2Class(find(ID2Class(:,2)==c),1);
    for j = jobsInClass(:)'
        jobsInClassArvClass = jobArvClass{1+j};
        jobsInClassRespT = jobRespT{1+j};
        if isempty(classResT{c})
            classResT{c} = jobsInClassRespT(jobsInClassArvClass==c);
        else
            classResT{c} = [classResT{c}; jobsInClassRespT(jobsInClassArvClass==c)];
        end
    end
end
end

