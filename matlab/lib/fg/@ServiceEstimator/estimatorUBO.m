function [estVal,fObjFun] = estimatorUBO(self, nodes)
% UBO Utilization-based optimization
% This demand estimator is based on the method proposed in:
%
% Liu, Z., Wynter, L., Xia, C. H. and Zhang, F.
% Parameter inference of queueing models for IT systems using end-to-end measurements
% Performance Evaluation, Elsevier, 2006.
%
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

%%
% rescale utilization to be mean number of busy servers
sn = self.model.getStruct;
V = zeros(sn.nclasses, size(nodes,2));
if ~iscell(nodes)
    nodes = {nodes};
end

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

        for chain=1:size(sn.visits, 1)
            V(r,n) = V(r, n) + sn.visits{chain}(node.index, r); 
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

[estVal, fObjFun] = ubo_data(avgU, avgR, avgA, V, self.options.iter_max, self.options.variant);
end

function [demEst,fObjFun] = ubo_data(cpuUtil, rAvgTimes, avgArvR, visits, ITERMAX, TYPE)
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


%% optimization program
N = size(cpuUtil,1); % number of experiments= size(cpuUtil,1); % number of experiments
L = nnz(~visits);

theta0 = zeros(M, R);

x = theta0';
x = x(:);

Aeq = diag(x);
Aeq(Aeq == -1) = 0;
Aeq(Aeq ~= 0) = 1;
beq = x;
beq(beq == -1) = 0;

XN = {}; 
HN = {};

As = visits;
As(As ~= 0) = 1;
As = diag(1 - As);
As(~any(As,2), :) = [];

for n=1:N
    nthAvgArvR = avgArvR(n, :);
    nthCpuUtil = cpuUtil(n, :);
    nthRespAvgTimes = rAvgTimes(n, :);

    % weight matrix
    w = nthAvgArvR./sum(nthAvgArvR);

    % H matrix
    Hs = zeros(M*R);
    Heps = 2 * M * ones(M);
    Hdelta = 2 * diag(w);
    nH = blkdiag(Hs, Heps, Hdelta);
    
    % A matrix

    arvV = diag(nthAvgArvR) * visits;
    cellArvV = num2cell(arvV', 1);
    cellArvV = cellfun(@diag, cellArvV, 'UniformOutput', false);
    Ae = cell2mat(cellArvV);

    cellAd = num2cell(visits, 2);
    cellAd = cellfun(@(x) x .* beta(n), cellAd, 'UniformOutput', false);
    Ad = blkdiag(cellAd{:});

    nA = [As, zeros(L, M), zeros(L, R);
        Ae, (-M*ones(M)), zeros(M, R);
        Ad, zeros(R, M), (-M*ones(R))];

    % b matrix
    nb = [zeros(L,1), nthCpuUtil, nthRespAvgTimes]';
    
    if TYPE == "BUNDLE"
        H(n, :, :) = nH;
        A(n, :, :) = nA;
        b(n, :) = nb;
    else
        if isempty(XN)
            [x, fObjFun]=fmincon(@nestedInitialObj, zeros(M * R + M + R, 1), [],[], nA, nb, zeros(M * R + M + R, 1), [],[], options);
        else
            [x, fObjFun]=fmincon(@nestedObj, zeros(M * R + M + R, 1), [],[], nA, nb, zeros(M * R + M + R, 1), [], [], options);
        end
        if fObjFun <= 100
            HN{end + 1} = nH;
            XN{end + 1} = x;
        end
    end

end

if TYPE == "BUNDLE"
    A = reshape(A, size(A,1)*size(nA,1), size(nA, 2));
    b = reshape(b, size(b,1)*size(nb,1), size(nb, 2));
    [x, fObjFun]=fmincon(@bundleObj, zeros(M * R + M + R, 1), [], [], A, b, XLB, XUB, options);
end

demEst = (reshape(x(1:M*R), M, R));

    function f = nestedInitialObj(x)
         f = 0.5 * x' * nH * x;
    end

    function f = nestedObj(x)
         y1 = norm(Aeq * x(1:M*R) - beq);
         xs = repmat({x}, 1, length(XN));
         cell = cellfun(@(x, xn, hn) hn * (x - xn), xs, XN, HN, 'UniformOutput', false);
         y2 = cell2mat(cell');

         f =  (0.5 * x' * nH * x) + y1 + sum(y2 .^ 2);
    end

    function f = bundleObj(x)
        f = 0.5 * x' * squeeze(sum(H, [1])) * x;
    end
    
end
