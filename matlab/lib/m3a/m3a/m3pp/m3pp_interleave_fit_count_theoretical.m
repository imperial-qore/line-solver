function [fit,m3pps] = m3pp_interleave_fit_count_theoretical(MMAP, t, tinf, mapping)
% Interleaves k M3PP to fit a multi-class trace with m classes.
% INPUT
% - MMAP: model to fit
% - t: finite time scale
% - tinf: near-infinite time scale
% - mapping: m x k binary matrix mapping the m classes to k m3pp (optional,
%            by default k = m)

% number of classes
m = size(MMAP,2)-2;

% by default, map one class per m3pp
%if nargin < 4
%    mapping = eye(m);
%end

% by default, quadratic
if nargin < 4
    delta = 0.75;
    SIGMA = mmap_count_mcov(MMAP,tinf);
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
%mapping

% check mapping
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
a = map_count_mean(MMAP,1);

% per-class rates
ac = mmap_count_mean(MMAP, 1);

% filters
f = cell(k,1);
for j = 1:k
    f{j} = mapping(:,j) == 1;
end

% total rates of each m3pp
av = zeros(k,1);
for j = 1:k
    av(j) = sum(ac(mapping(:,j) == 1));
end

% per-class rates within each m3pp
acc = cell(k,1);
for j = 1:k
    acc{j} = ac(f{j});
end

% per-class IDC(t) e IDC(inf)
vtc = mmap_count_var(MMAP, t);
vinfc = mmap_count_var(MMAP, tinf);
btc = vtc ./ (ac * t);
binfc = vinfc ./ (ac * tinf);

% compute covariance between classes
stc = mmap_count_mcov(MMAP, t);

% compute IDC(t) e IDC(inf) for each m3pp
btv = zeros(k,1);
binfv = zeros(k,1);
for j = 1:k
    % group classes that belong to j and classes that don't belong to j
    n = size(MMAP{1},1);
    mmap2 = cell(1,4);
    mmap2{1} = MMAP{1};
    mmap2{2} = MMAP{2};
    mmap2{3} = zeros(n,n);
    for i = 1:m
        if mapping(i,j)
            mmap2{3} = mmap2{3} + MMAP{2+i};
        end
    end
    mmap2{4} = mmap2{2} - mmap2{3};
    % IDC
    vtv = mmap_count_var(mmap2,t);
    vinfv = mmap_count_var(mmap2,tinf);
    btv(j) = vtv(1) ./ (av(j)*t);
    binfv(j) = vinfv(1) ./ (av(j)*tinf);
end

vtc = mmap_count_var(MMAP,t);
gtc = cell(k,1);
for j = 1:k
    % number of classes in this partition
    mj = sum(mapping(:,j));
    % compute
    n = 1;
    for i = 1:m
        if mapping(i,j) == 1
            % simple case, this is the only class
            if mj == 1
                gtc{j} = vtc(i);
                break;
            end
            % complex case
            m3pp3 = cell(1,5);
            m3pp3{1} = MMAP{1};
            m3pp3{2} = MMAP{2};
            % this class
            m3pp3{3} = MMAP{2+i};
            % other classes in the same partition
            m3pp3{4} = zeros(size(MMAP{1},1));
            for h = 1:m
                if h ~= i && mapping(h,j) == 1
                    m3pp3{4} = m3pp3{4} + MMAP{2+h};
                end
            end
            % classes in other partitions
            m3pp3{5} = m3pp3{2} - m3pp3{3} - m3pp3{4};
            % marginal covariance within the partition
            covji = mmap_count_mcov(m3pp3, t);
            gtc{j}(n) = vtc(i) + covji(1,2);
            % next
            n = n + 1;
        end
    end
end

[fit,m3pps] = m3pp_interleave_fit_count(av, btv, binfv, ...
    acc, gtc, t, tinf, mapping);

end