%{
%{
 % @file pfqn_ab_amva.m
 % @brief Akyildiz-Bolch AMVA method for multi-server BCMP networks.
%}
%}

%{
%{
 % @brief Akyildiz-Bolch AMVA method for multi-server BCMP networks.
 % @fn pfqn_ab_amva(D, N, V, nservers, sched, fcfsSchmidt, marginalProbMethod)
 % @param D Service time matrix (M x K).
 % @param N Population vector (1 x K).
 % @param V Visit ratio matrix (M x K).
 % @param nservers Number of servers at each station (M x 1).
 % @param sched Scheduling strategies for each station.
 % @param fcfsSchmidt Whether to use Schmidt formula for FCFS stations.
 % @param marginalProbMethod Method for marginal probability calculation ("ab" or "scat").
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return RN Residence times.
 % @return CN Cycle times.
 % @return XN System throughput.
 % @return totiter Total number of iterations.
%}
%}
function [QN,UN,RN,CN,XN,totiter] = pfqn_ab_amva(D,N,V,nservers,sched,fcfsSchmidt,marginalProbMethod)
% [QN,UN,RN,CN,XN,totiter] = PFQN_AB_AMVA(D,N,V,NSERVERS,SCHED,FCFSSCHMIDT,MARGINALPROBMETHOD)
%
% Akyildiz-Bolch AMVA method for multi-server BCMP networks.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 6
    fcfsSchmidt = false;
end
if nargin < 7
    marginalProbMethod = 'ab';
end

M = size(D, 1);
K = size(D, 2);

[QN, UN, RN, CN, XN, totiter] = ab_linearizer(K, M, N, nservers, sched, V, D, fcfsSchmidt, marginalProbMethod);

end

function [QN, UN, RN, CN, XN, totiter] = ab_linearizer(K, M, population, nservers, type, v, s, fcfsSchmidt, marginalProbMethod)
% Akyildiz-Bolch linearizer method for multi-server BCMP networks.

% Initialize queue length matrix
L = zeros(M, K);
for i = 1:M
    for r = 1:K
        L(i, r) = population(r) / M;
    end
end

% Initialize L_without_r matrices (cell array of M x K matrices)
lWithoutR = cell(1, K);
for r = 1:K
    lWithoutR{r} = zeros(M, K);
    for i = 1:M
        for t = 1:K
            if r == t
                lWithoutR{r}(i, t) = (population(r) - 1) / M;
            else
                lWithoutR{r}(i, t) = L(i, r);
            end
        end
    end
end

% Fractional changes matrix (initialized to 0)
D = zeros(M, K, K);

% STEP 1: Apply core at full population
[LUpdated, ~, ~, ~, ~, ~] = pfqn_ab_core(K, M, population, nservers, type, v, s, 100, D, L, fcfsSchmidt, marginalProbMethod);

% STEP 2: Apply core at N-e_k populations
for r = 1:K
    populationWithoutC = population;
    populationWithoutC(r) = population(r) - 1;

    lWithoutC = zeros(M, K);
    for j = 1:M
        for c = 1:K
            lWithoutC(j, c) = lWithoutR{c}(j, c);
        end
    end

    [retQ, ~, ~, ~, ~, ~] = pfqn_ab_core(K, M, populationWithoutC, nservers, type, v, s, 100, D, lWithoutC, fcfsSchmidt, marginalProbMethod);

    for j = 1:M
        for c = 1:K
            lWithoutR{c}(j, r) = retQ(j, c);
        end
    end
end

% STEP 3: Compute estimates of F_mk(N) and F_mk(N-e_j)
for i = 1:M
    for r = 1:K
        F_ir = LUpdated(i, r) / population(r);
        for t = 1:K
            if r == t
                divisor = population(r) - 1;
            else
                divisor = population(r);
            end
            if divisor ~= 0
                F_irt = lWithoutR{r}(i, t) / divisor;
            else
                F_irt = 0;
            end
            if isnan(F_irt)
                F_irt = 0;
            end
            D(i, r, t) = F_irt - F_ir;
        end
    end
end

% STEP 4: Apply core at full population using L values from step 1 and D values from step 3
[QN, UN, RN, CN, XN, totiter] = pfqn_ab_core(K, M, population, nservers, type, v, s, 100, D, LUpdated, fcfsSchmidt, marginalProbMethod);

end

function [QN, UN, W, CN, XN, totiter] = pfqn_ab_core(K, M, population, nservers, type, v, s, maxiter, D, lIn, fcfsSchmidt, marginalProbMethod)
% Akyildiz-Bolch core method for multi-server BCMP networks.

L = lIn;
tol = 1 / (4000 + 16 * sum(population));

totiter = 0;
W = zeros(M, K);

while totiter < maxiter
    F = zeros(M, K);
    lWithoutJ = zeros(M, K, K);

    for i = 1:M
        for r = 1:K
            if population(r) > 0
                F(i, r) = L(i, r) / population(r);
            else
                F(i, r) = 0;
            end
        end
    end

    for i = 1:M
        for r = 1:K
            for t = 1:K
                if r == t
                    scalar = population(r) - 1;
                else
                    scalar = population(r);
                end
                lWithoutJ(i, r, t) = scalar * (F(i, r) + D(i, r, t));
            end
        end
    end

    for i = 1:M
        for r = 1:K
            if type(i) == SchedStrategy.INF
                W(i, r) = s(i, r);
            elseif nservers(i) == 1
                totalQueueLengthKvecC = 0;
                for c = 1:K
                    totalQueueLengthKvecC = totalQueueLengthKvecC + lWithoutJ(i, c, r);
                end
                W(i, r) = s(i, r) * (1 + totalQueueLengthKvecC);
            elseif fcfsSchmidt && type(i) == SchedStrategy.FCFS
                waitTime = 0;
                nvec = pprod(population);
                while all(nvec >= 0)
                    if nvec(r) > 0
                        Bcn = getBcnForAB(s, i, r, nvec, K, nservers(i));
                        prob = getMarginalProb(oner(nvec, r), oner(population, r), population, lIn(i, r), K, i);
                        waitTime = waitTime + (Bcn * prob);
                    end
                    nvec = pprod(nvec, population);
                end
                if waitTime <= 1e-3
                    waitTime = 0;
                end
                W(i, r) = waitTime;
            else
                queueLength = 0;
                for j = 1:K
                    queueLength = queueLength + lWithoutJ(i, j, r);
                end
                numServers = floor(nservers(i));
                multiServerStationWeightedQueueLength = 0;
                if numServers > 1
                    populationWithoutR = population;
                    populationWithoutR(r) = populationWithoutR(r) - 1;
                    marginalProbs = findMarginalProbs(queueLength, numServers, populationWithoutR, r, marginalProbMethod);

                    for j = 1:(numServers-1)
                        if isfield(marginalProbs, num2str(j)) || isKey(marginalProbs, j)
                            prob_j = marginalProbs(j);
                        else
                            prob_j = 0;
                        end
                        multiServerStationWeightedQueueLength = multiServerStationWeightedQueueLength + prob_j * (numServers - j);
                    end
                end
                waitTime = (s(i, r) / numServers) * (1 + queueLength + multiServerStationWeightedQueueLength);
                W(i, r) = waitTime;
            end
        end
    end

    % Calculate cycle time for each class
    CN = zeros(1, K);
    for r = 1:K
        cycleTime = 0;
        for i = 1:M
            cycleTime = cycleTime + v(i, r) * W(i, r);
        end
        CN(r) = cycleTime;
    end

    % Calculate queue length L_ir = N_r * W_ir / C_r
    iterationQueueLength = zeros(M, K);
    for i = 1:M
        for r = 1:K
            if CN(r) > 0
                queueLength = population(r) * (v(i, r) * W(i, r) / CN(r));
            else
                queueLength = 0;
            end
            iterationQueueLength(i, r) = queueLength;
        end
    end

    % Check convergence
    maxDifference = 0;
    for i = 1:M
        for r = 1:K
            if population(r) > 0
                difference = abs(L(i, r) - iterationQueueLength(i, r)) / population(r);
                maxDifference = max(maxDifference, difference);
            end
        end
    end

    totiter = totiter + 1;
    L = iterationQueueLength;

    if maxDifference < tol
        break
    end
end

% Calculate throughput
XN = zeros(1, K);
for r = 1:K
    if W(1, r) > 0
        XN(r) = L(1, r) / W(1, r);
    else
        XN(r) = 0;
    end
end

% Calculate utilization
UN = zeros(M, K);
for i = 1:M
    for r = 1:K
        if s(i, r) > 0
            if type(i) == SchedStrategy.INF
                UN(i, r) = XN(r) * s(i, r);
            else
                UN(i, r) = (XN(r) * s(i, r)) / nservers(i);
            end
        else
            UN(i, r) = 0;
        end
    end
end

QN = L;

end

function Bcn = getBcnForAB(D, i, c, nvec, K, ns)
% Calculate Bcn for AB algorithm
Bcn = D(i, c);
if sum(nvec) > 1
    eps_val = 1e-12;
    sumVal = 0;
    for t = 1:K
        sumVal = sumVal + nvec(t) * D(i, t);
    end
    Bcn = Bcn + (max(0, sum(nvec) - ns) / max(ns * (sum(nvec) - 1), eps_val) * (sumVal - D(i, c)));
end
end

function prob = getMarginalProb(n, k, K_pop, L_jr, R, j)
% Compute marginal probability using binomial distribution
prob = 1;
for r = 1:R
    frac = L_jr / K_pop(r);
    if frac ~= 0 && K_pop(r) > 0 && n(r) >= 0 && n(r) <= K_pop(r)
        term1 = nchoosek(K_pop(r), n(r));
        term2 = frac^n(r);
        term3 = (1 - frac)^(K_pop(r) - n(r));
        prob = prob * (term1 * term2 * term3);
    end
end
end

function marginalProbs = findMarginalProbs(avgJobs, numServers, population, classIdx, marginalProbMethod)
% Finds marginal probabilities using specified method.

marginalProbs = containers.Map('KeyType', 'double', 'ValueType', 'double');

if strcmp(marginalProbMethod, 'scat')
    floorVal = floor(avgJobs);
    ceilVal = floorVal + 1;
    marginalProbs(floorVal) = ceilVal - avgJobs;
    marginalProbs(ceilVal) = avgJobs - floorVal;
    return
end

% AB method
ALPHA = 45.0;
BETA = 0.7;

w = weightFun(population, ALPHA, BETA);

floorVal = floor(avgJobs);
ceiling = floorVal + 1;
maxVal = min((2 * floorVal) + 1, numServers - 2);

for j = 0:maxVal
    if j <= floorVal
        lDist = floorVal - j;
        lowerVal = floorVal - lDist;
        upperVal = ceiling + lDist;
        if lDist > 25
            prob = 0;
        else
            if floorVal < population(classIdx)
                if upperVal ~= lowerVal
                    prob = w(floorVal+1, lDist+1) * ((upperVal - avgJobs) / (upperVal - lowerVal));
                else
                    prob = 0;
                end
            else
                prob = 0;
            end
        end
        marginalProbs(j) = prob;
    else
        uDist = j - ceiling;
        if uDist > 25
            marginalProbs(j) = 0;
        elseif j > population(classIdx) - 1 && uDist < 25
            if isKey(marginalProbs, population(classIdx) - 1)
                existingProb = marginalProbs(population(classIdx) - 1);
            else
                existingProb = 0;
            end
            if isKey(marginalProbs, floorVal - uDist)
                mp_floor_udist = marginalProbs(floorVal - uDist);
            else
                mp_floor_udist = 0;
            end
            newProb = existingProb + (w(floorVal+1, uDist+1) - mp_floor_udist);
            marginalProbs(population(classIdx) - 1) = newProb;
        else
            if isKey(marginalProbs, floorVal - uDist)
                mp_floor_udist = marginalProbs(floorVal - uDist);
            else
                mp_floor_udist = 0;
            end
            newProb = w(floorVal+1, uDist+1) - mp_floor_udist;
            marginalProbs(j) = newProb;
        end
    end
end

end

function w = weightFun(population, alpha, beta)
% Computes weight function for marginal probability calculation.

maxClassPopulation = max(population);

% Calculate scaling function PR
scalingFun = zeros(1, maxClassPopulation + 1);
if maxClassPopulation >= 1
    scalingFun(2) = alpha;  % scalingFun[1] in 0-indexed = scalingFun(2) in 1-indexed
    for n = 2:maxClassPopulation
        scalingFun(n+1) = beta * scalingFun(n);
    end
end

% Calculate weight function W (using 1-indexed matrices)
w = zeros(maxClassPopulation + 1, maxClassPopulation + 1);
w(1, 1) = 1.0;  % weightFun[0][0] = 1.0

for l = 1:maxClassPopulation
    for j = 0:(l-1)
        w(l+1, j+1) = w(l, j+1) - (w(l, j+1) * scalingFun(l+1)) / 100.0;
    end
    sumVal = 0;
    for j = 0:(l-1)
        sumVal = sumVal + w(l+1, j+1);
    end
    w(l+1, l+1) = 1 - sumVal;
end

end
