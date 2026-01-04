%{
%{
 % @file pfqn_schmidt_ext.m
 % @brief Extended Schmidt MVA algorithm with queue-aware alpha corrections.
%}
%}

%{
%{
 % @brief Extended Schmidt MVA algorithm with queue-aware alpha corrections.
 % @fn pfqn_schmidt_ext(D, N, S, sched)
 % @param D Service demand matrix (M x R).
 % @param N Population vector (1 x R).
 % @param S Number of servers per station (M x 1).
 % @param sched Scheduling discipline per station.
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return CN Cycle times.
 % @return T Results table.
%}
%}
function [XN,QN,UN,CN,T] = pfqn_schmidt_ext(D,N,S,sched)
% [XN,QN,UN,CN,T] = PFQN_SCHMIDT_EXT(D,N,S,SCHED)
%
% Extended Schmidt MVA algorithm with queue-aware alpha corrections.
% A queue-aware version of the Schmidt algorithm that precomputes alpha values
% for improved accuracy in networks with class-dependent FCFS scheduling.
%
% Reference: R. Schmidt, "An approximate MVA algorithm for exponential,
% class-dependent multiple server stations," Performance Evaluation,
% vol. 29, no. 4, pp. 245-254, 1997.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[M, R] = size(D);
closedClasses = 1:R;
XN = zeros(1, R);
UN = zeros(M, R);
CN = zeros(M, R);
QN = zeros(M, R);

C = length(closedClasses);
Dc = D(:, closedClasses);
Nc = N(closedClasses);

% Precompute alphas for FCFS stations
alphas = cell(M, R);
for i = 1:M
    if sched(i) == SchedStrategy.FCFS
        % Check if class-dependent
        classIndependent = true;
        for r = 2:R
            if D(i, r) ~= D(i, 1)
                classIndependent = false;
                break
            end
        end

        if ~classIndependent
            for r = 1:R
                % Create modified problem with tagged customer
                D_mod = zeros(M, R + 1);
                N_mod = zeros(1, R + 1);

                for k = 1:R
                    if k == r
                        N_mod(k) = N(k) - 1;
                    else
                        N_mod(k) = N(k);
                    end
                    for j = 1:M
                        D_mod(j, k) = D(j, k);
                        if i == j
                            D_mod(j, R + 1) = D(j, r);
                        else
                            D_mod(j, R + 1) = 0;
                        end
                    end
                end
                N_mod(R + 1) = 1;

                % Create extended sched array
                sched_mod = [sched(:)', SchedStrategy.FCFS];

                % Run basic Schmidt to get alpha values
                [~, ~, U_alpha, ~, ~] = pfqn_schmidt(D_mod, N_mod, S, sched_mod);
                alphas{i, r} = U_alpha;
            end
        end
    end
end

% Compute products for hashing
prods = zeros(1, C);
for r = 1:C
    prods(r) = prod(Nc(1:r-1) + 1);
end

% Start at nc=(0,...,0)
kvec = pprod(Nc);

% Initialize L and Pc
L = cell(1, M);
Pc = cell(1, M);
for ist = 1:M
    % Always initialize L for ALL stations (queue lengths needed everywhere)
    L{ist} = zeros(R, prod(1 + Nc));
    switch sched(ist)
        case SchedStrategy.INF
            % No Pc needed for infinite server
        case SchedStrategy.PS
            if ~all(S(ist, :) == 1)
                % Multi-server PS needs state probabilities
                Pc{ist} = zeros(1 + sum(Nc), prod(1 + Nc));
            end
        case SchedStrategy.FCFS
            classIndependent = all(D(ist, :) == D(ist, 1));
            if classIndependent
                if ~all(S(ist, :) == 1)
                    Pc{ist} = zeros(1 + sum(Nc), prod(1 + Nc));
                end
            else
                % Class-dependent needs vector probabilities
                Pc{ist} = zeros(prod(1 + Nc), prod(1 + Nc));
            end
    end
end

x = zeros(C, prod(1 + Nc));
w = zeros(M, C, prod(1 + Nc));

for ist = 1:M
    if ~isempty(Pc{ist})
        Pc{ist}(1 + 0, hashpop(kvec, Nc, C, prods)) = 1.0;
    end
end

% Population recursion
while all(kvec >= 0) && all(kvec <= Nc)
    kprods = zeros(1, C);
    for r = 1:C
        kprods(r) = prod(kvec(1:r-1) + 1);
    end

    for ist = 1:M
        for c = 1:C
            hkvec = hashpop(kvec, Nc, C, prods);
            hkvec_c = hashpop(oner(kvec, c), Nc, C, prods);

            if length(S(ist, :)) == 1
                ns = S(ist);
            else
                ns = S(ist, c);
            end

            if kvec(c) > 0
                switch sched(ist)
                    case SchedStrategy.INF
                        w(ist, c, hkvec) = D(ist, c);

                    case SchedStrategy.PS
                        if ns == 1
                            % Sum queue lengths of ALL classes (arrival theorem)
                            totalQueueLength = sum(L{ist}(:, hkvec_c));
                            w(ist, c, hkvec) = Dc(ist, c) * (1 + totalQueueLength);
                        else
                            % Sum queue lengths of ALL classes (arrival theorem)
                            totalQueueLength = sum(L{ist}(:, hkvec_c));
                            w(ist, c, hkvec) = (Dc(ist, c) / ns) * (1 + totalQueueLength);
                            for j = 1:ns-1
                                % Pc{ist}(j,...) = Pr(j-1 jobs) due to MATLAB 1-indexing
                                w(ist, c, hkvec) = w(ist, c, hkvec) + (ns - 1 - (j - 1)) * Pc{ist}(j, hkvec_c) * (Dc(ist, c) / ns);
                            end
                        end

                    case SchedStrategy.FCFS
                        classIndependent = all(D(ist, :) == D(ist, 1));

                        if ~classIndependent
                            if ns == 1
                                % Sum queue lengths of ALL classes (arrival theorem)
                                totalQueueLength = sum(L{ist}(:, hkvec_c));
                                w(ist, c, hkvec) = Dc(ist, c) * (1 + totalQueueLength);
                            else
                                nvec = pprod(kvec);
                                while all(nvec >= 0)
                                    if nvec(c) > 0
                                        % Use Nc and prods for hash (matches Java)
                                        hnvec_c = hashpop(oner(nvec, c), Nc, C, prods);

                                        % Use extended Bcn with alpha
                                        if Nc(c) > 1 && ~isempty(alphas{ist, c})
                                            Bcn = getBcnExt(alphas{ist, c}, D, ist, c, nvec, C, ns, R);
                                        else
                                            Bcn = getBcn(D, ist, c, nvec, C, ns);
                                        end

                                        w(ist, c, hkvec) = w(ist, c, hkvec) + Bcn * Pc{ist}(hnvec_c, hkvec_c);
                                    end
                                    nvec = pprod(nvec, kvec);
                                end
                            end
                        else
                            % Class-independent: treat like PS
                            if ns == 1
                                % Sum queue lengths of ALL classes (arrival theorem)
                                totalQueueLength = sum(L{ist}(:, hkvec_c));
                                w(ist, c, hkvec) = Dc(ist, c) * (1 + totalQueueLength);
                            else
                                % Sum queue lengths of ALL classes (arrival theorem)
                                totalQueueLength = sum(L{ist}(:, hkvec_c));
                                w(ist, c, hkvec) = (Dc(ist, c) / ns) * (1 + totalQueueLength);
                                for j = 1:ns-1
                                    % Pc{ist}(j,...) = Pr(j-1 jobs) due to MATLAB 1-indexing
                                    w(ist, c, hkvec) = w(ist, c, hkvec) + (ns - 1 - (j - 1)) * Pc{ist}(j, hkvec_c) * (Dc(ist, c) / ns);
                                end
                            end
                        end
                end
            end
        end
    end

    % Compute throughputs
    for c = 1:C
        sumW = sum(w(1:M, c, hkvec));
        if sumW > 0
            x(c, hkvec) = kvec(c) / sumW;
        else
            x(c, hkvec) = 0;
        end
    end
    x(isnan(x)) = 0;

    % Update queue lengths
    for ist = 1:M
        for c = 1:C
            L{ist}(c, hkvec) = x(c, hkvec) * w(ist, c, hkvec);
        end

        if length(S(ist, :)) == 1
            ns = S(ist);
        else
            ns = S(ist, c);
        end

        switch sched(ist)
            case SchedStrategy.PS
                if ns > 1
                    for n = 1:min(S(ist), sum(kvec))
                        for c = 1:C
                            if kvec(c) > 0
                                hkvec_c = hashpop(oner(kvec, c), Nc, C, prods);
                                Pc{ist}(1 + n, hkvec) = Pc{ist}(1 + n, hkvec) + Dc(ist, c) * (1/n) * x(c, hkvec) * Pc{ist}(1 + (n - 1), hkvec_c);
                            end
                        end
                        Pc{ist}(1 + 0, hkvec) = max(eps, 1 - sum(Pc{ist}(1 + (1:min(S(ist), sum(kvec))), hkvec)));
                    end
                end

            case SchedStrategy.FCFS
                classIndependent = all(D(ist, :) == D(ist, 1));

                if ~classIndependent
                    nvec = pprod(kvec);
                    nvec = pprod(nvec, kvec);
                    sumOfAllProbs = 0;
                    while all(nvec >= 0)
                        % Use Nc and prods for hash (matches Java)
                        hnvec = hashpop(nvec, Nc, C, prods);

                        prob = 0;
                        for r = 1:C
                            if nvec(r) > 0
                                % Use Nc and prods for hash (matches Java)
                                hnvec_c = hashpop(oner(nvec, r), Nc, C, prods);
                                hkvec_c = hashpop(oner(kvec, r), Nc, C, prods);

                                if Nc(r) > 1 && ~isempty(alphas{ist, r})
                                    Bcn = getBcnExt(alphas{ist, r}, D, ist, r, nvec, C, ns, R);
                                else
                                    Bcn = getBcn(D, ist, r, nvec, C, ns);
                                end

                                capacity_inv = 1 / sum(nvec);
                                x_ir = x(r, hkvec);
                                prob_c = Pc{ist}(hnvec_c, hkvec_c);
                                classProb = Bcn * capacity_inv * x_ir * prob_c;
                                prob = prob + classProb;
                            end
                        end
                        Pc{ist}(hnvec, hkvec) = prob;
                        sumOfAllProbs = sumOfAllProbs + prob;
                        nvec = pprod(nvec, kvec);
                    end
                    Pc{ist}(1 + 0, hkvec) = max(1e-12, 1 - sumOfAllProbs);
                else
                    % Class-independent multi-server
                    if ns > 1
                        K_j = 0;
                        for r = 1:R
                            % Assuming all classes visit if D > 0
                            if D(ist, r) > 0
                                K_j = K_j + Nc(r);
                            end
                        end

                        meanQueueLength = 0;
                        for r = 1:R
                            meanQueueLength = meanQueueLength + L{ist}(r, hkvec);
                        end

                        for n = 1:(min(ns, sum(kvec))-1)
                            if K_j > 0 && n <= K_j
                                frac = meanQueueLength / K_j;
                                x1 = nchoosek(K_j, n);
                                x2 = frac^n;
                                x3 = (1 - frac)^(K_j - n);
                                prob = x1 * x2 * x3;
                                Pc{ist}(1 + n, hkvec) = prob;
                            end
                        end

                        sum1 = 0;
                        sum2 = 0;
                        for r = 1:R
                            sum1 = sum1 + D(ist, r) * x(r, hkvec);
                        end
                        for n = 0:(ns-1)
                            sum2 = sum2 + (ns - n) * Pc{ist}(1 + n, hkvec);
                        end
                        term = (sum1 + sum2) / ns;
                        Pc{ist}(1 + 0, hkvec) = max(1e-12, 1 - term);
                    end
                end
        end
    end

    kvec = pprod(kvec, Nc);
end

% Extract final results
hkvecFinal = hkvec;

% Throughput
XN(closedClasses) = x(1:C, hkvecFinal);
if M > 1
    XN = repmat(XN, M, 1);
end

% Utilization
for m = 1:M
    for c = 1:C
        UN(m, c) = (D(m, c) * XN(c)) / S(m);
    end
end

% Response time
CN(1:M, closedClasses) = w(1:M, 1:C, hkvecFinal);

% Queue length
for ist = 1:M
    QN(ist, closedClasses) = L{ist}(closedClasses, hkvecFinal);
end

T = table(XN, CN, QN, UN);

end

function Bcn = getBcn(D, i, c, nvec, C, ns)
% Standard Bcn calculation
Bcn = D(i, c);
if sum(nvec) > 1
    eps_val = 1e-12;
    sumVal = 0;
    for t = 1:C
        sumVal = sumVal + nvec(t) * D(i, t);
    end
    Bcn = Bcn + (max(0, sum(nvec) - ns) / max(ns * (sum(nvec) - 1), eps_val) * (sumVal - D(i, c)));
end
end

function Bcn = getBcnExt(u, D, i, c, nvec, C, ns, R)
% Extended Bcn with alpha corrections
weightedProb = 0;

% u is the utilization matrix from the alpha computation (M x R+1)
% The last column (R+1) is the tagged customer utilization
totalNonPinnedTime = sum(u(i, 1:C)) - u(i, C + 1);

if totalNonPinnedTime > 0
    for s = 1:C
        if D(i, s) > 0
            prob = u(i, s) / totalNonPinnedTime;
            weightedProb = weightedProb + (prob / D(i, s));
        end
    end
end

if weightedProb > 0
    meanInterdepartureTime = 1.0 / (ns * weightedProb);
else
    meanInterdepartureTime = 0;
end

Bcn = D(i, c);

if sum(nvec) > 1
    Bcn = Bcn + (max(0, sum(nvec) - ns) * meanInterdepartureTime);
end

if isnan(Bcn) || isinf(Bcn)
    Bcn = 0;
end
end
