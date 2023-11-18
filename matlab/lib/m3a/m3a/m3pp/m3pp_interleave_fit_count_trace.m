function [fit,m3pps] = m3pp_interleave_fit_count_trace(T, A, t, tinf, mapping)
% Interleaves k M3PP to fit a multi-class trace with m classes.
% INPUT
% - T: inter-arrival time
% - A: labels
% - t: finite time scale
% - tinf: near-infinite time scale
% - mapping: m x k binary matrix mapping the m classes to k m3pp (optional,
%            by default k = m)


% number of classes
L = unique(A);
m = length(L);

% by default choose time scales as follows
if nargin < 3
    t = 10*mean(T);
    tinf = max(10*t, (sum(T)-T(1)) / 100);
end
% compute counting process at resolution t1 and tinf
Nt = mtrace_iat2counts(T, A, t);
Ninf = mtrace_iat2counts(T, A, tinf);

% by default, %%map one class per m3pp
if nargin < 5
    %    mapping = eye(m);
    delta = 0.75;
    SIGMA = cov(Ninf);
    pool = 1:m;
    J=1;
    variance = diag(SIGMA);
    while length(pool)>0
        k = pool(maxpos(variance(pool)));
        mapping(k,J)=1;
        pool(find(k == pool))=[];
        for h=pool
            if SIGMA(k,h)/sqrt(SIGMA(k,k)*SIGMA(h,h)) >= delta
                mapping(h,J)=1;
                pool(find(h == pool))=[];
            end
        end
        J = J+1;
    end
    
end
% check mapping consistency
if size(mapping,1) ~= m
    error('Number of classes does not match mapping.');
end

% number of m3pps
k = size(mapping,2);

% check mapping
rowsum = sum(mapping,2);
if max(rowsum) > 1 || min(rowsum) < 1
    error('Invalid mapping');
end

fprintf('Fitting %d classes with %d M3PP(2,m_j) processes\n', m, k);

% rate
a = 1/mean(T);

% per-class rates
pc = zeros(m,1);
for j = 1:m
    pc(j) = sum(A==L(j))/length(A);
end
ac = pc * a;

% filters
f = cell(k,1);
for j = 1:k
    f{j} = mapping(:,j) == 1;
end

% total rates of each m3pp
av = zeros(k,1);
for j = 1:k
    av(j) = sum(ac(f{j}));
end

% per-class rates within each m3pp
acc = cell(k,1);
for j = 1:k
    acc{j} = ac(f{j});
end

% compute IDC(t) e IDC(inf) for each m3pp
btv = zeros(k,1);
binfv = zeros(k,1);
for j = 1:k
    % counting processes for the j-th m3pp
    Ntj = zeros(size(Nt,1),1);
    Ninfj = zeros(size(Ninf,1),1);
    for i = 1:m
        % sum classes fitted by this m3pp
        if mapping(i,j)
            Ntj = Ntj + Nt(:,i);
            Ninfj = Ninfj + Ninf(:,i);
        end
    end
    % IDC
    btv(j) = var(Ntj)/(av(j)*t);
    binfv(j) = var(Ninfj)/(av(j)*tinf);
end

% compute per-class variance plus marginal covariance within each m3pp
gtc = cell(k,1);
for j = 1:k
    % number of classes fitted with this m3pp
    mj = sum(mapping(:,j));
    % filter classes
    fj = mapping(:,j) == 1;
    % per-class variance of this m3pp
    vtj = var(Nt(:,fj))';
    % per class variance + marginal covariance of this m3pp
    stj = zeros(mj,1);
    h = 1;
    for i = 1:m
        % if class i is mapped to the j-th m3pp
        if mapping(i,j) == 1
            covtji = cov(Nt(:,i), sum(Nt(:,fj),2)-Nt(:,i));
            stj(h) = covtji(1,2);
            h = h + 1;
        end
    end
    % vector of per-class variance plus margina covariance
    gtc{j} = vtj + stj;
end

[fit,m3pps] = m3pp_interleave_fit_count(av, btv, binfv, ...
    acc, gtc, t, tinf, mapping);

end