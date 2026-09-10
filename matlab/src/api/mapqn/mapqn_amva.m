function [X, Qq, U, pik, ES] = mapqn_amva(mu, D0s, D1s, N)
% [X, QQ, U, PIK, ES] = MAPQN_AMVA(MU, D0S, D1S, N)
%
% Horizontal-cut mean value analysis of a closed multiclass network made of
% an exponential infinite-server station (think rate MU(r) for class r) and
% one FCFS single-server station whose class-r service process is the MAP
% (D0S{r}, D1S{r}). The MAP of class r moves only while a class-r job is in
% service and is frozen otherwise, the convention of SolverCTMC.
%
% The recursion walks the population lattice n <= N in lexicographic order and
% at every point solves ONE linear R x R system. Its unknowns are the per-phase
% means Q_r^k = E[n_r 1{k}] over the joint phase k = (k_1..k_R), the busy laws
% U_r^k = P[serving r, k], the phase law pi_k and the throughputs X_r. Exact
% relations: the joint phase balance, the class marginals U_r = X_r E[S_r]
% theta_r, and the per-class horizontal cut (generator balance of n_r 1{k}) of
% Casale-Smirni, DSN 2009. Closures: the product busy law theta_r(k_r) times
% the post-completion laws of the frozen MAPs (which solves the phase balance
% identically), the service-age closure of the cross term E[n_r 1{serving s}]
% (class r accumulates at its throughput over the elapsed class-s service,
% whose mean given the phase is theta_s(-D0_s)^{-1} / theta_s), and Little's
% law resolved by arrival phase with the exact FCFS response of the queue
% composition seen at population n - e_r (the multiclass arrival theorem).
% K_r = 1 for every class reproduces multiclass FCFS MVA on class means.
%
% Inputs:
%   mu   - 1 x R think rates (exponential infinite server)
%   D0s  - 1 x R cell, D0s{r} the K_r x K_r hidden-transition matrix of class r
%   D1s  - 1 x R cell, D1s{r} the K_r x K_r completion matrix of class r
%   N    - 1 x R populations (a class with N(r) = 0 is absent)
%
% Outputs:
%   X    - 1 x R class throughputs
%   Qq   - 1 x R mean queue lengths at the MAP station (in service included)
%   U    - 1 x R busy probabilities of the server per class, X(r) E[S_r]
%   pik  - 1 x prod(K_r) joint phase law at population N (class R fastest)
%   ES   - 1 x R mean service times
%
% Reference:
%   G. Casale, E. Smirni, "MAP-AMVA: Approximate Mean Value Analysis of Bursty
%   Systems", IEEE/IFIP DSN 2009, pp. 409-418 (the horizontal cut).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

R = numel(N);
mu = mu(:)';
N = round(N(:)');
Ks = zeros(1, R);
for r = 1:R
    Ks(r) = size(D0s{r}, 1);
end
K = prod(Ks);
stride = zeros(1, R);
for r = 1:R
    stride(r) = prod(Ks(r+1:end));           % class R is the fastest index
end
krof = zeros(K, R);                           % 1-based phase of class r in joint phase k
for k = 1:K
    for r = 1:R
        krof(k, r) = mod(floor((k - 1) / stride(r)), Ks(r)) + 1;
    end
end
Z = 1 ./ mu;
Ntot = sum(N);

G = cell(1, R); th = cell(1, R); phi = cell(1, R); T = cell(1, R);
age = cell(1, R); Ainv = cell(1, R); ES = zeros(1, R); abar = zeros(1, R);
for r = 1:R
    Kr = Ks(r);
    G{r} = D0s{r} + D1s{r};
    th{r} = stationary_(G{r});
    ES(r) = 1 / (th{r} * sum(D1s{r}, 2));
    phi{r} = th{r} * D1s{r} * ES(r);          % post-completion phase law
    s = (-D0s{r}) \ ones(Kr, 1);              % mean time to the next completion from phase k
    P = (-D0s{r}) \ D1s{r};                   % embedded phase transition at completions
    Tr = zeros(Ntot + 2, Kr);                 % Tr(j+1,k) = e_k' (I + P + .. + P^(j-1)) s
    acc = zeros(Kr, 1); v = s;
    for j = 1:Ntot + 1
        acc = acc + v; Tr(j + 1, :) = acc'; v = P * v;
    end
    T{r} = Tr;
    w = th{r} / (-D0s{r});                    % theta (-D0)^{-1}
    age{r} = w ./ th{r};                      % mean elapsed service given the phase
    abar(r) = sum(w);
    Ainv{r} = inv(G{r} - mu(r) * eye(Kr));
end

% joint-phase shapes: idle law F and class-r busy law u(r,:)
F = ones(1, K); u = ones(R, K);
for k = 1:K
    for r = 1:R
        F(k) = F(k) * phi{r}(krof(k, r));
        for s = 1:R
            if s == r
                u(r, k) = u(r, k) * th{s}(krof(k, s));
            else
                u(r, k) = u(r, k) * phi{s}(krof(k, s));
            end
        end
    end
end

% population lattice in lexicographic order: n - e_r always precedes n
lstride = zeros(1, R);
for r = 1:R
    lstride(r) = prod(N(r+1:end) + 1);
end
L = prod(N + 1);
Qs = zeros(L, R, K); pis = zeros(L, K); Xs = zeros(L, R);
pis(1, :) = F;
for l0 = 1:L - 1
    l = l0 + 1;
    n = zeros(1, R);
    for r = 1:R
        n(r) = mod(floor(l0 / lstride(r)), N(r) + 1);
    end
    present = find(n >= 1);
    % arrival-theorem conditionals at n - e_r, class responses, age closure
    b = cell(1, R); bN = cell(1, R); Rk = zeros(R, K);
    for r = present
        lp = l0 - lstride(r) + 1;
        pip = pis(lp, :); Qp = reshape(Qs(lp, :, :), R, K);
        br = zeros(R, K); pos = pip > 0;
        br(:, pos) = Qp(:, pos) ./ repmat(pip(pos), R, 1);
        b{r} = br;
        for k = 1:K
            acc = 0;
            for t = 1:R
                acc = acc + tat_(T{t}, br(t, k) + double(t == r), krof(k, t));
            end
            Rk(r, k) = acc;
        end
        bN{r} = zeros(R, K);
        for t = 1:R
            if t == r || n(t) == 0
                continue
            end
            Xt = Xs(lp, t); Qt = sum(Qp(t, :));
            if Xt > 0
                W = max(Qt / Xt - abar(r), 0);
                for k = 1:K
                    bN{r}(t, k) = Xt * min(W + age{r}(krof(k, r)), n(t) / Xt);
                end
            end
        end
    end
    % the cut, linear in X: Q_r = c0{r} + sum_s X(s) c1{r,s}
    c0 = cell(1, R); c1 = cell(R, R);
    for r = present
        c0{r} = axis_(-mu(r) * n(r) * F, Ainv{r}, r, krof, stride, Ks(r), K);
        for s = 1:R
            term = -mu(r) * n(r) * ES(s) * (u(s, :) - F);
            if s == r
                term = term + ES(r) * axis_(u(r, :), D1s{r}, r, krof, stride, Ks(r), K);
            elseif n(s) >= 1
                W = ES(s) * u(s, :) .* bN{s}(r, :);
                term = term + axis_(W, G{r}, r, krof, stride, Ks(r), K) ...
                    - axis_(W, G{s}, s, krof, stride, Ks(s), K);
            end
            c1{r, s} = axis_(term, Ainv{r}, r, krof, stride, Ks(r), K);
        end
    end
    % Little's law by arrival phase: one R x R solve
    Mm = eye(R); v = zeros(R, 1);
    for r = present
        v(r) = n(r) - mu(r) * sum((n(r) * F - c0{r}) .* Rk(r, :));
        Mm(r, r) = Z(r);
        for s = present
            Mm(r, s) = Mm(r, s) + mu(r) * sum((n(r) * ES(s) * (u(s, :) - F) - c1{r, s}) .* Rk(r, :));
        end
    end
    X = (Mm \ v)';
    pik = F;
    for s = 1:R
        pik = pik + X(s) * ES(s) * (u(s, :) - F);
    end
    pik = max(pik, 0); pik = pik / sum(pik);
    for r = present
        Ur = X(r) * ES(r) * u(r, :);
        Qr = c0{r};
        for s = 1:R
            Qr = Qr + X(s) * c1{r, s};
        end
        % project onto Q >= U keeping the flow-balance total n_r - X_r/mu_r
        Wr = max(Qr - Ur, 0); tot = max(n(r) - X(r) / mu(r) - sum(Ur), 0);
        if sum(Wr) > 0
            Qr = Ur + Wr * (tot / sum(Wr));
        else
            Qr = Ur;
        end
        Qs(l, r, :) = Qr;
    end
    pis(l, :) = pik; Xs(l, :) = X;
end
X = Xs(L, :);
Qq = sum(reshape(Qs(L, :, :), R, K), 2)';
U = X .* ES;
pik = pis(L, :);
end

function th = stationary_(G)
K = size(G, 1);
A = [G'; ones(1, K)]; bvec = [zeros(K, 1); 1];
th = (A \ bvec)';
th = max(th, 0); th = th / sum(th);
end

function out = axis_(V, M, r, krof, stride, Kr, K)
% out(k) = sum_h V(k with k_r -> h) M(h, k_r): M applied along the class-r axis
out = zeros(1, K);
for k = 1:K
    kr = krof(k, r);
    base = (k - 1) - (kr - 1) * stride(r);
    acc = 0;
    for h = 1:Kr
        acc = acc + V(base + (h - 1) * stride(r) + 1) * M(h, kr);
    end
    out(k) = acc;
end
end

function val = tat_(Tr, b, kr)
% linear interpolation of the response table at a fractional number of jobs ahead
j0 = floor(b); j0 = max(0, min(j0, size(Tr, 1) - 2)); f = min(max(b - j0, 0), 1);
val = (1 - f) * Tr(j0 + 1, kr) + f * Tr(j0 + 2, kr);
end
