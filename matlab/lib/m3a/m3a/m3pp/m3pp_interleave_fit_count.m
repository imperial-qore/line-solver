function [fit,m3pps] = m3pp_interleave_fit_count(av, btv, binfv, ...
                                                 acc, gtcc, t, tinf, ...
                                                 mapping, reorder)
% Fits k second-order M3PP[m_j] and interleaves them into a
% M3PP[m] of order k+1, with m = \sum_j=1^k m_j.
%
% INPUT
%  av:      vector of length k with the per-process rate
%  btv:     vector of length k with the per-process IDC(t)
%  binfv:   vector of length k with the per-process IDC(inf)
%  acc:     cell-array of length k; j-th element is a vector of length m_j
%           with the per-class rates in the j-th M3PP
%  gtcc:    cell-array of length k; j-th element is a vector of length m_j
%           with the per-class variance + marginal_covariance in the j-th
%           M3PP
%  t:       finite time scale
%  tinf:    near-infinite time scale
%  mapping: binary m x k matrix mapping classes to processes
%  reorder: (optional, default = 1), enables reordering
% 
% OUTPUT
%  fit:     result of the interleaving, M3PP[m] of order k+1
%  m3pps:   cell-array with the fitted and interleaved second-order
%           M3PP[m_j]

if nargin < 9
    reorder = 1;
end

% get number of partitions
k = length(av);

% number of classes for each partition
mv = zeros(k,1);
for j = 1:k
    mv(j) = length(acc{j});
end

% total number of classes
m = sum(mv);

% check consistency
for j = 1:k
    if abs(av(j)-sum(acc{j})) > 1e-8
        error('Inconsistent per-class rates in the %-th partition', j); 
    end
end

% vector of lower bounds for the upper off-diagonal elements
uv = zeros(k,1);
% vector of upper bounds for the upper off-diagonal elements
dv = zeros(k,1);

% compute bounds for the upper off-diagonal elment of each mmpp(2)
for j = 1:k    
    if btv(j) >= binfv(j)
        d = compute_d(btv(j), binfv(j)*(1+1e-4), t);
    else
        d = compute_d(btv(j), binfv(j), t);
    end
    z = (binfv(j)-1)*d^3*av(j);
    u = d*z/(2*av(j)^2*d^2+z);
    uv(j) = u;
    dv(j) = d;
end

% find feasible values for off diagonal elements
if ~reorder
    [r,tran] = compute_feasible_interleave(uv,dv);
    order = 1:k;
else
    [r,tran,order] = compute_feasible_interleave_reorder(uv,dv);
end

% compute position of each class
cpos = zeros(m,1);
for c = 1:m
    % process modeling this class
    j = find(mapping(c,:) == 1);
    % position of the process
    jpos = find(order == j);
    % compute class position
    cpos(c) = 0;
    for hpos = 1:(jpos-1)
        h = order(hpos);
        cpos(c) = cpos(c) + sum(mapping(:,h));
    end
    cpos(c) = cpos(c) + sum(mapping(1:c,j));
end

% fit m3pp[m_j] processes
m3pps = cell(k,1);
for o = 1:k
    % get process in o-th position
    j = order(o);
    % fit underlying mmpp(2)
    %tran
    if ~tran(j)
        r1 = r(1,j);
        r2 = r(2,j);
    else
        r1 = r(2,j);
        r2 = r(1,j);
    end
    d = r1 + r2;
    z = (binfv(j)-1)*d^3*av(j);
    delta = sqrt(z/(2*r1*r2));
    l2 = av(j) - r2/d * delta;
    l1 = l2 + delta;
    if tran(j)
        tmp = r1;
        r1 = r2;
        r2 = tmp;
        tmp = l1;
        l1 = l2;
        l2 = tmp;
    end
    mmpp = cell(1,2);
    mmpp{1}(1,2) = r1;
    mmpp{1}(2,1) = r2;
    mmpp{2}(1,1) = l1;
    mmpp{2}(2,2) = l2;
    drates = -sum(mmpp{1}+mmpp{2},2);
    for h = 1:2
        mmpp{1}(h,h) = drates(h);
    end
    % get vectors of per-class characteristics
    acj = acc{j};
    gtcj = gtcc{j};
    % fit m3pp(2,mj)
    m3pps{o} = m3pp2m_fit_count_approx_ag_multiclass(mmpp, acj, gtcj, t);
end

% total number of classes
m = 0;
for j = 1:k
    m = m + length(acc{j});
end

% compute interleaving
fit_unordered = m3pp2m_interleave(m3pps);

% re-order output classes
fit = cell(1,2+m);
fit{1} = fit_unordered{1};
fit{2} = fit_unordered{2};
for i = 1:m
    fit{2+i} = fit_unordered{2+cpos(i)};
end

% check feasibility and normalize
if ~mmap_isfeasible(fit, 1e-10)
    warning('Fitted MMAP infeasibile at tolerance 1e-10\n');
end
fit = mmap_normalize(fit);

% find mapping of classes to the m3pps
J = zeros(m,1);
O = zeros(m,1);
for i = 1:m
    j = find(mapping(i,:) == 1);
    o = sum(mapping(1:i,j));
    J(i) = j;
    O(i) = o;
end

% total rate
a = sum(av);
fa = map_count_mean(fit,1);

% per class rates
ac = zeros(m,1);
for i = 1:m
    ac(i) = acc{J(i)}(O(i));
end
fac = mmap_count_mean(fit,1);

% compare characteristics
fprintf('Rate: input = %f, %f\n', a, fa);
for i = 1:m
    fprintf('Class %d, rate: input = %f, output = %f\n', i, ac(i), fac(i));
end
fbtv = zeros(k,1);
for j = 1:k
    m3pp2 = cell(1,4);
    m3pp2{1} = fit{1};
    m3pp2{2} = fit{2};
    m3pp2{3} = zeros(k+1,k+1);
    for i = 1:m
        if mapping(i,j) == 1
            m3pp2{3} = m3pp2{3} + fit{2+i};
        end
    end
    m3pp2{4} = m3pp2{2}-m3pp2{3};
    fvj = mmap_count_var(m3pp2,t);
    fbtv(j) = fvj(1)/(av(j)*t);
end
for j = 1:k
    fprintf('Partition %d, IDC(%.2f): input = %f, output = %f\n', j, t, btv(j), fbtv(j));
end
fbinfv = zeros(k,1);
fbinfv_limit = zeros(k,1);
for j = 1:k
    m3pp2 = cell(1,4);
    m3pp2{1} = fit{1};
    m3pp2{2} = fit{2};
    m3pp2{3} = zeros(k+1,k+1);
    for i = 1:m
        if mapping(i,j) == 1
            m3pp2{3} = m3pp2{3} + fit{2+i};
        end
    end
    m3pp2{4} = m3pp2{2}-m3pp2{3};
    fvj = mmap_count_var(m3pp2,tinf);
    fbinfv(j) = fvj(1)/(av(j)*tinf);
    fvj = mmap_count_var(m3pp2, 1e6);
    fbinfv_limit(j) = fvj(1)/(av(j)*1e6);
end
for j = 1:k
    fprintf('Partition %d, IDC(inf = %.2f): input = %f, output = %f, output_limit = %f\n', j, tinf, binfv(j), fbinfv(j), fbinfv_limit(j));
end
fvtc = mmap_count_var(fit,t);
for j = 1:k
    n = 1;
    for i = 1:m
        if mapping(i,j) == 1
            m3pp3 = cell(1,5);
            m3pp3{1} = fit{1};
            m3pp3{2} = fit{2};
            m3pp3{3} = fit{2+i};
            m3pp3{4} = zeros(k+1,k+1);
            for h = 1:m
                if h ~= i && mapping(h,j) == 1
                    m3pp3{4} = m3pp3{4} + fit{2+h};
                end
            end
            m3pp3{5} = m3pp3{2} - m3pp3{3} - m3pp3{4};
            covji = mmap_count_mcov(m3pp3, t);
            fprintf('gamma(%d): ', i);
            fprintf('input = %f, output = %f\n', gtcc{j}(n), fvtc(i) + covji(1,2));
            n = n + 1;
        end
    end
end

    % solve non-linear equation to fit underlying mmpp(2)
    function d = compute_d(bt1,binf,t1)
        if ~(binf>bt1 && bt1>1)
            error('No solution, infeasible IDC(t): IDC(%.2f) = %.3f, IDC(inf) = %.3f\n',...
                   t1, bt1, binf);
        end
        % d = r1 + r2
        g = (binf-1)/(binf-bt1);
        % according to WolframAlpha
        % the solution can be written as (ProductLog[-c e^(-c)]+c)/t1
        z = -g*exp(-g);
        w = fsolve(@(w) z-w*exp(w),  1, optimset('Display','none',...
                                                 'MaxIter',10000,...
                                                 'MaxFunEvals',10000,...
                                                 'TolFun', 1.0e-12,...
                                                 'TolX',1.0e-12));
        d = (w + g)/t1;
    end

end