%{
%{
 % @file pfqn_pam.m
 % @brief Hsieh-Lam Proportional Approximation Methods (PAMB/PAMI/PAMT).
%}
%}

function [XN,QN,UN,RN]=pfqn_pam(L,N,Z,variant)
%{
%{
 % @brief Hsieh-Lam Proportional Approximation Methods (PAMB/PAMI/PAMT).
 %
 % C. T. Hsieh, S. S. Lam, "PAM - A noniterative approximate solution method
 % for closed multichain queueing networks", ACM SIGMETRICS Perform. Eval.
 % Rev. 16(1), 1988. The three variants are NONITERATIVE: the queue lengths
 % are seeded by the proportion of a class demand that falls at each centre,
 %
 %   E_ck = D_ck / sum_i D_ci,   Q_ck(N) = E_ck N_c,
 %
 % and the MVA equations are then unrolled a fixed number of times.
 %
 %   'pamb'  seed, then the last MVA step (eqs. 2.29-2.34)
 %   'pami'  PAMB, then throughputs scaled down wherever a centre would be
 %           driven past full utilization (eqs. 2.35 and step 4)
 %   'pamt'  seed at N - 1_i - 1_j, then the last TWO MVA steps, then the
 %           PAMI utilization capping (eqs. 2.36-2.40)
 %
 % The seed spreads the whole class population over the queueing centres and
 % ignores Z, exactly as published: PAM buys speed, not accuracy.
 %
 % @fn pfqn_pam(L, N, Z, variant)
 % @param L Service demand matrix (stations x classes).
 % @param N Population vector.
 % @param Z Think time vector.
 % @param variant 'pamb' (default), 'pami' or 'pamt'.
 % @return XN System throughput.
 % @return QN Mean queue lengths.
 % @return UN Utilization.
 % @return RN Residence times.
%}
%}

[M,R]=size(L);
if nargin<3 || isempty(Z)
    Z = zeros(1,R);
end
Z = sum(Z,1);
if nargin<4 || isempty(variant)
    variant = 'pamb';
end
variant = lower(variant);

% E_ck, the share of the class-r demand served at station ist (eq. 2.29)
E = zeros(M,R);
Ltot = sum(L,1);
for r=1:R
    if Ltot(r) > 0
        E(:,r) = L(:,r)/Ltot(r);
    end
end
Q = E .* repmat(N,M,1);   % Q_ck(N) = E_ck N_c (eq. 2.30)

RN = zeros(M,R);
XN = zeros(1,R);
switch variant
    case 'pamt'
        for i=1:R
            % Q_ck(N - 1_i - 1_j) = Q_ck(N) - E_ck [(c==i) + (c==j)] (eq. 2.36)
            Qmi = zeros(M,R);   % Q_jk(N - 1_i), filled by the inner sweep
            for j=1:R
                Qij = Q - E .* repmat((1:R)==i,M,1) - E .* repmat((1:R)==j,M,1);
                Rj = L(:,j) .* (1 + sum(Qij,2));            % eq. (2.37)
                beta = double(i==j);
                nj = N(j) - beta;
                if nj > 0
                    Xj = nj/(sum(Rj) + Z(j));               % eqs. (2.38)-(2.39)
                else
                    Xj = 0;
                end
                Qmi(:,j) = Xj*Rj;                           % eq. (2.40)
            end
            RN(:,i) = L(:,i) .* (1 + sum(Qmi,2));           % step 3(b)
            if N(i) > 0
                XN(i) = N(i)/(sum(RN(:,i)) + Z(i));
            end
        end
    otherwise % 'pamb' and the first step of 'pami'
        for r=1:R
            % Q_jk(N - 1_r) = Q_jk(N) - E_jk [j == r] (eq. 2.31)
            Qmr = Q - E .* repmat((1:R)==r,M,1);
            RN(:,r) = L(:,r) .* (1 + sum(Qmr,2));           % eq. (2.32)
            if N(r) > 0
                XN(r) = N(r)/(sum(RN(:,r)) + Z(r));         % eqs. (2.33)-(2.34)
            end
        end
end

if strcmp(variant,'pami') || strcmp(variant,'pamt')
    % scale a class down when it would drive a centre it visits past U = 1
    U = sum(L .* repmat(XN,M,1), 2);                        % eq. (2.35)
    for r=1:R
        visited = L(:,r) ~= 0;
        if any(visited)
            S = max(U(visited));
            if S > 1
                XN(r) = XN(r)/S;
            end
        end
    end
end

QN = repmat(XN,M,1) .* RN;
UN = repmat(XN,M,1) .* L;
end
