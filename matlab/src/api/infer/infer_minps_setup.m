function demandEst = infer_minps_setup(data, initSample, sampleSize, V, model, node)
% MAIN_MINPS setups the input data for the MINPS estimation method and calls it
%
% Copyright (c) 2012-2014, Imperial College London 
% All rights reserved.


R = size(data,2) - 1;
sampleNumber = zeros(1,R+1);
for k = 1:R
    sampleNumber(k) = size(data{3,k},1);
end

% remove classes without samples
newR = sum(sampleNumber>0);
data2 = cell(6,newR+1);
r=1;
sampleNumber(R+1) = 1;
for k=1:R+1
    if sampleNumber(k)>0;
        for j=1:6
            data2{j,r} = data{j,k};
        end
        r=r+1;
    end
end
data = data2;
R = newR;
% get queue length at arrival times
qls = infer_get_qlen_arrival(data);

rt = [];    %response times
class = []; % job classes
ql = [];
at = [];
for k = 1:R
    rt = [rt; data{4,k}];
    class = [class; k*ones(size(data{4,k},1), 1)];
    ql = [ql; qls{k} ];
    at = [at; data{3,k}/1000];
end

% sort data to include data from all classes in the data set
allTimes = [at rt class ql];
allTimes = sortrows(allTimes,1);

at = allTimes(:,1);
rt = allTimes(:,2);
class = allTimes(:,3);
ql = allTimes(:,4:end);

if sampleSize == 0
    sampleSize = size(ql,1);
end

% select sample set
firstSample = initSample;
finalSample = initSample+sampleSize-1;
sampleSet = firstSample:finalSample;

qlExp = ql(sampleSet,:);
rtExp = rt(sampleSet);
classExp = class(sampleSet);
numClassExp = accumarray(classExp(:), 1, [R 1])';

% remove samples with zero response times
valid = rtExp > 0;
if any(~valid)
    rtExp = rtExp(valid);
    classExp = classExp(valid);
    qlExp = qlExp(valid,:);
end

% Estimate number of threads 
Wexp = max(sum(ql,2));
% Estimate mean number of busy processors
numNotProcExp = Wexp - sum(sum(ql,1),2)/size(ql,1); 
numNotProcExp = numNotProcExp/R;
% Estimate think-time rates 
lambda = zeros(1,R);
for k = 1:R
    lambda(k) = min(1E6, (numClassExp(k)/(at(finalSample) + rt(finalSample) - at(firstSample)) )/numNotProcExp);  
    
    if lambda(k) < 0 
        lambda(k) = 1E6;
    end
end

demandEst = infer_minps(model, node, rtExp, classExp, qlExp);
