function [estVal,fObjFun] = estimatorUBO2(self, node)
% UBO Utilization-based optimization
% This demand estimator is based on the method proposed in:
%
% Daniel A. Menasc
% Computing Missing Service Demand Parameters for Performance Models
%
% Copyright (c) 2012-2022, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

% This estimator at present can handle estimation on a single resource.

%%
% rescale utilization to be mean number of busy servers
sn = self.model.getStruct;

if isfinite(node.getNumberOfServers())
    U = self.getAggrUtil(node);
    if ~isempty(U)
        avgU = U.data * node.getNumberOfServers();
    end
end

% obtain per class metrics
for r=1:sn.nclasses
    avgArvR{r} = self.getArvR(node, self.model.classes{r});
    if isempty(avgArvR{r})
        error('Arrival rate data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
    else
        avgArvR{r} = avgArvR{r}.data;
    end
    avgRespT{r} = self.getRespT(node, self.model.classes{r});
    if isempty(avgRespT{r})
        error('Response time data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
    else
        avgRespT{r} = avgRespT{r}.data;
    end
    
end

try
    avgA = cell2mat(avgArvR);
    avgR = cell2mat(avgRespT);
    sumUr = sum(cell2mat(avgUtil),2);
catch me
    switch me.identifier
        case 'MATLAB:catenate:dimensionMismatch'
            error('Sampled metrics have different number of samples, use interpolate() before starting this estimation algorithm.');
    end
end

[estVal, fObjFun] = ubo_data(avgU, avgR, avgA, self.options.iter_max);
estVal = estVal(:)';
end

% ubo procedure based on the comon data format
function [demEst,fObjFun] = ubo_data(cpuUtil, rAvgTimes, avgArvR, ITERMAX)
a = isnan(cpuUtil);
if sum(a) > 0
    disp('NaN values found for CPU Utilization. Removing NaN values.');
    cpuUtil = cpuUtil(a == 0);
    rAvgTimes = rAvgTimes(a == 0,:);
    avgArvR = avgArvR(a == 0,:);
end

a = sum(avgArvR,2) == 0;
if sum(a) > 0
    disp('Removing sampling intervals with zero throughput for all request types.');
    cpuUtil = cpuUtil(a == 0);
    rAvgTimes = rAvgTimes(a == 0,:);
    avgArvR = avgArvR(a == 0,:);
end


%% number of resources
M = 1;
%% number of classes
R = size(rAvgTimes,2);

beta = repmat(1./(1-cpuUtil),1,R);

%% initial point
% x0(r) is the mean service demand of class r (visits are assumed unitary)
x0 = rand(M,R).*max(rAvgTimes); % randomize service demand in [0,max(avgRTime)] for each class
%% options
options = optimset();
options.Display = 'off';
options.LargeScale = 'off';
options.MaxIter =  ITERMAX;
options.MaxFunEvals = 1e10;
options.MaxSQPIter = 5000;
options.TolCon = 1e-8;
options.Algorithm = 'interior-point';

XLB = x0*0 + options.TolCon; % lower bounds on x variables
XUB = repmat(max(rAvgTimes), M,1); % upper bounds on x variables

T0 = tic; % needed for outfun

%% optimization program
N = size(cpuUtil,1); % number of experiments= size(cpuUtil,1); % number of experiments

[demEst, fObjFun]=fmincon(@objfun,x0,[],[],[],[],XLB,XUB,[],options);

    function f = objfun(x)
        d = repmat(x,N,1);
        f = 0;
        for c=1:R
            rC = 0;
            for i=1:M
                denomR = 0;
                for cN=1:R
                    denomR = denomR + d(1,cN,i)*avgArvR(1, cN); 
                end
                rC = rC + d(1,c,i)/(1-denomR);
            end
            f = f + (rC-rAvgTimes(1,c)).^2;
        end
    end

end
