function pfqn_sens_validate()
% PFQN_SENS_VALIDATE  Numerically validate the queue-length moment / demand-
% derivative identities for closed product-form (BCMP) queueing networks.
%
% Notation (single-server load-independent stations, plus optional delay Z):
%   Q_{i,r}(N)          mean queue length of class r at station i, pop. N
%                       -> pfqn_mva(L,N,Z), QN(i,r)
%   Q_{i,s}^{+k}(N)     queue length of class s at ORIGINAL station i in the
%                       network obtained by adding one replica of station k
%                       (a station with identical demands L(k,:))
%                       -> pfqn_mva([L;L(k,:)],N,Z), read row i
%   D_{i,r} = L(i,r)    service demand
%   dQ_{k,r}/dD_{i,s}   -> pfqn_sens(L,N,Z), dQ(k,r,pL(i,s))
%   N - 1_r             population with one class-r job removed
%
% The three "moment" identities and the three demand-derivative formulas from
% the prompt are each checked over all admissible index tuples (i,k,r,s) on a
% set of random closed networks. Analytic derivatives (pfqn_sens) are also
% cross-checked against central finite differences of pfqn_mva.

rng(12345);

% ---- test models: {L, N, Z} -------------------------------------------------
models = {};
models{end+1} = struct('L',[0.6 0.3; 0.4 0.7; 0.5 0.2],'N',[3 2],'Z',[0 0]);
models{end+1} = struct('L',[1.0 0.5; 0.3 0.9; 0.7 0.6],'N',[2 3],'Z',[0.4 0.8]);
models{end+1} = struct('L',rand(4,3)+0.2,             'N',[2 2 2],'Z',[0 0 0]);
models{end+1} = struct('L',rand(3,2)+0.1,             'N',[4 3],  'Z',[0.5 0]);

tolIdent = 1e-8;   % pure MVA identities: near machine precision
tolDeriv = 1e-6;   % vs analytic sens (both from MVA recursion)
tolFD    = 1e-4;   % analytic sens vs finite difference

maxErr = zeros(1,7);   % [ident1 ident2 ident3 der_diag der_offclass der_offstat sens_vs_FD]

for mm = 1:numel(models)
    L = models{mm}.L; N = models{mm}.N; Z = models{mm}.Z;
    [M,R] = size(L);

    [~,QN] = pfqn_mva(L,N,Z);                 % base Q_{i,r}(N)
    sens   = pfqn_sens(L,N,Z);                % analytic derivatives
    pL = @(i,s) (i-1)*R + s;                  % pfqn_sens L-parameter index

    % cache Q^{+k}(N-1_r) for every (k,r) with N_r>=1
    Qp = cell(M,R);                           % Qp{k,r} is M x R (rows=orig stat)
    for k = 1:M
        for r = 1:R
            if N(r) >= 1
                Qp{k,r} = Qplus(L,N,Z,k,r);
            end
        end
    end

    for r = 1:R
        if N(r) < 1, continue; end
        for s = 1:R
            for k = 1:M
                for i = 1:M
                    % ---------- moment identity 1 ----------
                    % Q_{i,s}^{+k}(N-1_r) Q_{k,r} = Q_{k,r}^{+i}(N-1_s) Q_{i,s}
                    if N(s) >= 1
                        lhs = Qp{k,r}(i,s) * QN(k,r);
                        rhs = Qp{i,s}(k,r) * QN(i,s);
                        maxErr(1) = max(maxErr(1), relerr(lhs,rhs));

                        % ---------- moment identity 2 ----------
                        % Q_{i,s}^{+k}(N-1_r) Q_{k,r} D_{k,s} D_{i,r}
                        %   = Q_{i,r}^{+k}(N-1_s) Q_{k,s} D_{i,s} D_{k,r}
                        % (prompt wrote the superscript as +i; the exact
                        %  identity requires +k, forced by identities 1 and 3)
                        lhs = Qp{k,r}(i,s)*QN(k,r)*L(k,s)*L(i,r);
                        rhs = Qp{k,s}(i,r)*QN(k,s)*L(i,s)*L(k,r);
                        maxErr(2) = max(maxErr(2), relerr(lhs,rhs));

                        % ---------- moment identity 3 ----------
                        % Q_{i,s}^{+k}(N-1_r) Q_{k,r} D_{k,s} D_{i,r}
                        %   = Q_{k,s}^{+i}(N-1_r) Q_{i,r} D_{i,s} D_{k,r}
                        lhs = Qp{k,r}(i,s)*QN(k,r)*L(k,s)*L(i,r);
                        rhs = Qp{i,r}(k,s)*QN(i,r)*L(i,s)*L(k,r);
                        maxErr(3) = max(maxErr(3), relerr(lhs,rhs));
                    end

                    % ---------- demand-derivative formulas ----------
                    % analytic LHS = D_{i,s} dQ_{k,r}/dD_{i,s}
                    dAna = L(i,s) * sens.dQ(k,r,pL(i,s));
                    % finite-difference cross-check of the analytic derivative
                    dFD  = L(i,s) * fd_dQ(L,N,Z,k,r,i,s);
                    maxErr(7) = max(maxErr(7), relerr(dAna,dFD));

                    if i==k && r==s
                        % D_{k,r} dQ_{k,r}/dD_{k,r}
                        %   = Q_{k,r}(1 + 2 Q_{k,r}^{+k}(N-1_r) - Q_{k,r})
                        rhs = QN(k,r)*(1 + 2*Qp{k,r}(k,r) - QN(k,r));
                        maxErr(4) = max(maxErr(4), relerr(dAna,rhs));
                    elseif i==k && r~=s
                        % D_{k,s} dQ_{k,r}/dD_{k,s}
                        %   = Q_{k,r}(2 Q_{k,s}^{+k}(N-1_r) - Q_{k,s})
                        rhs = QN(k,r)*(2*Qp{k,r}(k,s) - QN(k,s));
                        maxErr(5) = max(maxErr(5), relerr(dAna,rhs));
                    elseif i~=k
                        % D_{i,s} dQ_{k,r}/dD_{i,s}
                        %   = Q_{k,r}(Q_{i,s}^{+k}(N-1_r) - Q_{i,s})
                        rhs = QN(k,r)*(Qp{k,r}(i,s) - QN(i,s));
                        maxErr(6) = max(maxErr(6), relerr(dAna,rhs));
                    end
                end
            end
        end
    end
end

names = {'moment identity 1        ', ...
         'moment identity 2        ', ...
         'moment identity 3        ', ...
         'deriv i=k,r=s (diagonal) ', ...
         'deriv i=k,r~=s (offclass)', ...
         'deriv i~=k     (offstat) ', ...
         'sens vs finite-difference'};
tols  = [tolIdent tolIdent tolIdent tolDeriv tolDeriv tolDeriv tolFD];

fprintf('\n=== pfqn_sens formula validation (max relative error) ===\n');
allok = true;
for f = 1:7
    ok = maxErr(f) <= tols(f);
    allok = allok && ok;
    fprintf('  %s : %10.3e  (tol %7.1e)  %s\n', names{f}, maxErr(f), tols(f), passfail(ok));
end
fprintf('\n%s\n\n', ternary(allok,'ALL FORMULAS VALIDATED','SOME FORMULAS FAILED'));
if ~allok
    error('pfqn_sens_validate:mismatch','one or more formulas exceeded tolerance');
end
end

% -------------------------------------------------------------------------
function Q = Qplus(L,N,Z,k,r)
% Q_{i,s}^{+k}(N-1_r): queue lengths at the ORIGINAL M stations when a replica
% of station k is added, evaluated at population N with one class-r job removed.
Np = N; Np(r) = Np(r) - 1;
M  = size(L,1);
Lp = [L; L(k,:)];
[~,QNp] = pfqn_mva(Lp,Np,Z);
Q = QNp(1:M,:);
end

% -------------------------------------------------------------------------
function d = fd_dQ(L,N,Z,k,r,i,s)
% central finite difference of Q_{k,r}(N) w.r.t. L(i,s)
h = 1e-6 * max(1,abs(L(i,s)));
Lp = L; Lp(i,s) = Lp(i,s) + h;
Lm = L; Lm(i,s) = Lm(i,s) - h;
[~,Qp] = pfqn_mva(Lp,N,Z);
[~,Qm] = pfqn_mva(Lm,N,Z);
d = (Qp(k,r) - Qm(k,r)) / (2*h);
end

% -------------------------------------------------------------------------
function e = relerr(a,b)
e = abs(a-b) / max([1, abs(a), abs(b)]);
end

function s = passfail(ok), if ok, s='PASS'; else, s='FAIL'; end, end
function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
