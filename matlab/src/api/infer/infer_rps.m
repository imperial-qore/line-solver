function demandEst = infer_rps(rt, class, ql, V)
% INFER_RPS Regression for Processor Sharing (RPS) demand estimation
%
% RT:       response time samples (column vector with all samples)
% CLASS:    class of each request sample (column vector)
% QL:       queue length samples (matrix with R columns containing the
%           number of jobs of each class observed by each sample)
% V:        number of cores/servers
%
% Based on mean-value analysis for PS stations. For a single server:
%
%   E[R_r] = E[D_r] * E[Q_bar_A]
%
% where Q_bar_A is the total number of jobs seen upon admission,
% including the arriving job itself. For V > 1 processors, the
% queue length is split equally among the servers:
%
%   E[R_r] = E[D_r] * E[Q_bar_A] / V
%
% The demand E[D_r] is estimated via non-negative least squares
% regression of response times against Q_bar_A / V.
%
% Note: RPS assumes equal splitting among all V servers, which can
% lead to overestimation under low loads when fewer than V servers
% are actually busy. ERPS improves on this by using the mean number
% of busy servers instead of V.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

R = max(class);

demandEst = zeros(1, R);
for r = 1:R
    idx = (class == r);
    respTimes = rt(idx);
    % Q_bar_A: total jobs seen upon admission including the arriving job
    qBarA = sum(ql(idx, :), 2) + 1;
    demandEst(r) = lsqnonneg(qBarA / V, respTimes);
end
end
