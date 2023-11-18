function [estVal,fObjFun] = estimatorUBO3(self, nodes)
% UBO Utilization-based optimization
% This demand estimator is based on the method proposed in:
%
% Liu, Z., Wynter, L., Xia, C. H. and Zhang, F.
% Parameter inference of queueing models for IT systems using end-to-end measurements
% Performance Evaluation, Elsevier, 2006.
%
% Copyright (c) 2012-2022, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

%%
% rescale utilization to be mean number of busy servers
sn = self.model.getStruct;
for n=1:size(nodes,2)
    node = nodes{n};
    if isfinite(node.getNumberOfServers())
        U = self.getAggrUtil(node);
        if ~isempty(U)
            avgU{n} = U.data * node.getNumberOfServers();
        end
    end


    % obtain per class metrics
    for r=1:sn.nclasses
        avgArvR{n, r} = self.getArvR(node, self.model.classes{r});
        if isempty(avgArvR{n,r})
            error('Arrival rate data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
        else
            avgArvR{n,r} = avgArvR{n,r}.data;
        end
        avgRespT{n,r} = self.getRespT(node, self.model.classes{r});
        if isempty(avgRespT{n,r})
            error('Response time data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
        else
            avgRespT{n,r} = avgRespT{n,r}.data;
        end
        
    end
end


try
    avgA = cell2mat(avgArvR);
    avgA = reshape(avgA, size(avgA,1)/size(avgArvR, 1), size(avgArvR, 1), sn.nclasses);
    avgR = cell2mat(avgRespT);
    avgR = reshape(avgR, size(avgR,1)/size(avgRespT, 1), size(avgRespT, 1), sn.nclasses);
    avgU = cell2mat(avgU);
    sumUr = sum(avgU,2);
catch me
    switch me.identifier
        case 'MATLAB:catenate:dimensionMismatch'
            error('Sampled metrics have different number of samples, use interpolate() before starting this estimation algorithm.');
    end
end

[estVal, fObjFun] = ubo_data(avgU, avgR, avgA, self.options.iter_max);
end

% ubo procedure based on the comon data format
function [demEst,fObjFun] = ubo_data(cpuUtil, rAvgTimes, avgArvR, ITERMAX)
a = sum(isnan(cpuUtil), 2);
if sum(a) > 0
    disp('NaN values found for CPU Utilization. Removing NaN values.');
    cpuUtil = cpuUtil(a == 0, :);
    rAvgTimes = rAvgTimes(a == 0,:,:);
    avgArvR = avgArvR(a == 0,:,:);
end

a = sum(sum(avgArvR,3), 2) == 0;
if sum(a) > 0
    disp('Removing sampling intervals with zero throughput for all request types.');
    cpuUtil = cpuUtil(a == 0, :, :);
    rAvgTimes = rAvgTimes(a == 0, :, :);
    avgArvR = avgArvR(a == 0,:, :);
end

%% number of resources
M = size(cpuUtil, 2);
%% number of classes
R = size(rAvgTimes,3);

beta = repmat(1./(1-cpuUtil),1,1,R);

%% initial point
% x0(i, r) is the mean service demand of class r at station i(visits are assumed unitary)
x0 = rand(M,R).*squeeze(max(rAvgTimes,[],1)); % randomize service demand in [0,max(avgRTime)] for each class
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
XUB = squeeze(max(rAvgTimes)); % upper bounds on x variables

T0 = tic; % needed for outfun

%% optimization program
N = size(cpuUtil,1); % number of experiments= size(cpuUtil,1); % number of experiments
epsi = cpuUtil;
deltaj = cpuUtil;
w = avgArvR./repmat(sum(avgArvR,3), 1, 1, R);
[demEst, fObjFun]=fmincon(@objfun,x0,[],[],[],[],XLB,XUB,[],options);

    function f = objfun(x)
        d =  repmat(reshape(x, 1, size(x,1), size(x,2)), N, 1,1);
        epsi = sum(d.*avgArvR,3) - cpuUtil;
        deltaj = d.*beta - rAvgTimes;
        f = sum(sum(sum(w.*deltaj.^2, 3), 2)) + sum(sum(epsi.^2, 2));
    end
end
