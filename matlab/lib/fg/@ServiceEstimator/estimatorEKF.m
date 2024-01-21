function [estVal] = estimatorEKF(self, node)
% Extended Kalman Filter resource demand estimator 
% This demand estimator is based on the method proposed in:
%
% Kumar, Dinesh and Tantawi, Asser and Zhang, Li
% Real-Time Performance Modeling for Adaptive Software Systems 2009
%
% Copyright (c) 2012-2024, Imperial College London
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
    
    [estVal] = ekf_data(avgU, avgR, avgA, node.getNumberOfServers(), self.options.iter_max);
    estVal = estVal(:)';
end

% ekf procedure based on the comon data format
function [demEst] = ekf_data(cpuUtil, rAvgTimes, avgArvR, numServers, ITERMAX)
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
    
    %% initial state
    % x(r) is the mean service demand of class r (visits are assumed unitary)
    x = rand(R,1).*max(rAvgTimes); % randomize service demand in [0,max(avgRTime)] for each class
    p = diag(x.^2);

    % taken from Zheng, Woodside and Litoiu, Performance Model Estimation and Tracking Using Optimal Filters
    %pCovNoise = diag(0.0109, 0.0374, 0.0745, 0.0000154); 

    mCovNoise = diag(repmat(0.01, 1, R + 1)); 
    pCovNoise = eye(size(x,1)).*0.001;

    %% optimization program
    N = size(cpuUtil,1); % number of experiments= size(cpuUtil,1); % number of experiments
    stepBound = 0.6;
    a_min = zeros(size(x));
    a_max = zeros(size(x)) + inf;

    for n=1:min(N, ITERMAX)
        % predict
        % xn = Fn xn-1
        x_n = x;
        
        % Pn = Fn Pn-1 F'n + Qn
        P_n = p + pCovNoise;

        % update
        stepUtil = cpuUtil(n);
        stepArvR = avgArvR(n, :);
        stepResponse = rAvgTimes(n, :);

        % y~n =  zn - h(xn)
        H_n = getJacobian(x_n, R, stepUtil, stepArvR, numServers);
        z = getMeasurement(R, stepUtil, stepResponse);
        z_n = getPredictedMeasurement(x_n, R, stepArvR, numServers);
        y_n = z_n - z;

        % Sn =  Hn Pn H'n + Rn
        H_nT = H_n';
        S_n = H_n * P_n * H_nT + mCovNoise;

        % Kalman Gain Kn = Pn H'n S^-1n 
        K_n = P_n * H_nT * inv(S_n);

        % xnn = xn  + Kn y~n 
        x = x_n + K_n * y_n;
        xlower = (stepBound * a_min) + (1-stepBound)*x;
        xupper = (stepBound * a_max) + (1-stepBound)*x;
        x = min(xupper, max(xlower, x));
        if (sum(x) < 0)
            x = x * -1;
        end
        % Pnn = (I - Kn Hn) Pn
        p = (eye(size(x)) - (K_n * H_n))  * P_n;
        
    end
    [demEst] = x(1:M);
end

function [Hx] = getJacobian(x, C, util, arriv, P)
    Hx = zeros(C+1, C);
    for c=1:C
        Hx(c,c) = (1 - util + (arriv(c)*x(c))/P) / (1-util)^2;
        Hx(C+1,c) = arriv(c)/P;
    end
    [Hx] = Hx;
end

function [h] = getMeasurement(C, util, responseTimes)
    h = zeros(C+1,1);
    for c=1:C
        h(c) = responseTimes(c);
    end
    h(C+1) = util;
    [h] = h;
end

function [h] = getPredictedMeasurement(x, C, arriv, P)
    h = zeros(C+1,1);
    u = sum(x.*arriv)/P;
    for c=1:C
        h(c) = (x(c)/(1-u));
    end
    h(C+1) = u;
    [h] = h;
end
