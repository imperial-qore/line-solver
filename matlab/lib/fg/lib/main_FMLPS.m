function demandEst = main_FMLPS(data, initSample, sampleSize, V)
% MAIN_FMLPS setups the input data for the FMLPS estimation method and calls it
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
qls = getQLArrival(data);

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

% select sample set
firstSample = initSample;
finalSample = initSample+sampleSize-1;
sampleSet = firstSample:finalSample;

qlExp = ql(sampleSet,:);
rtExp = rt(sampleSet);
classExp = class(sampleSet);
numClassExp = hist(classExp, 1:R);

% remove samples with zero response times
rtzero = rtExp ==0;
if(sum(rtzero)>0)
    rtExp = rtExp(rtExp>0);
    classExp = classExp(rtExp>0);
    qlExp = qlExp(rtExp>0,:);
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

demandEst = des_FMLPS(rtExp, classExp, qlExp,lambda, V, Wexp); 
