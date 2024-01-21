function estVal = estimatorERPS(self, node)
% (rt,class,ql,nCores)
% DES_ERPS implements the ERPS demand estimation method
% RT:       response times samples.
%           column vector with all response time samples
% CLASS:    class of the request samples
%           column vector with the class of each sample
% QL:       queue length samples
%           matrix with R columns containing the number of jobs
%           of each class observed by each sample
% NCORES:   number of cores
%
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

%%
sn = self.model.getStruct;

if node.schedStrategy ~= SchedStrategy.PS
    error('The ERPS method is available only for processor sharing stations.');
end

nodeId = self.model.getNodeIndex(node);

% obtain per class metrics
R = sn.nclasses;
for r=1:R
    avgRespT{r} = self.getRespT(node, self.model.classes{r});        
    arvlEvent = Event(EventType.ID_ARV, node, self.model.classes{r}); %class-r arrival at node 
    avgAQLen{r} = self.getAggrQLen(node, arvlEvent); % aggregate queue-length
    
    if isempty(avgRespT{r})
        error('Response time data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
    else
        avgRespT{r} = avgRespT{r}.data;
    end
    if isempty(avgAQLen{r})
        error('Arrival queue-length data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
    else
        avgAQLen{r} = avgAQLen{r}.data;
    end
    if min(avgAQLen{r})<1 
        error('Arrival queue-length cannot be less than 1 as it must include the arriving job.');
    end
end

% estimate average number of busy cores
busyCores = sum(avgAQLen{1,1},2);
for r = 2:R
    busyCores = [busyCores;sum(avgAQLen{1,r},2)];
end
avgBusyCores = min(mean(busyCores), node.getNumberOfServers);

estVal = zeros(1,R);
% regression analysis to estimate mean demands
for r=1:R
    respTimes = avgRespT{1,r};
    totalQL = sum(avgAQLen{1,r},2)/avgBusyCores;
    estVal(r) = lsqnonneg(totalQL, respTimes)';
end
estVal = estVal(:)';
end