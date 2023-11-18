function [fit] = m3pp2m_fit_count_theoretical(mmap, method, t, tinf)
% Fits the theoretical characteristics of a MMAP(n,m) with a M3PP(2,m).
% INPUT:
% - mmap: the MMAP(n,m) to fit with a M3PP(2,m)
% - method: 'exact_delta', 'approx_delta', 'approx_cov' or 'approx_ag'
% - t: (optional) finite time scale
% - tinf: (optional) near-infinite time scale

if nargin == 1
    method = 'approx_delta';
end

m = size(mmap,2)-2;

if strcmp(method,'approx_cov') && m > 2
   error('Approximate covariance fitting only supported for two classes.');
end

if nargin < 3
    t1 = 1;
    t2 = 10;
    tinf = 1e4;
    t3 = 10;
else
    t1 = t;
    t2 = t;
    t3 = t;
end

% joint-process charactersitics
a = map_count_mean(mmap,t1)/t1;
bt1 = map_count_var(mmap,t1)/(a*t1);
bt2 = map_count_var(mmap,t2)/(a*t2);
binf = map_count_var(mmap,tinf)/(a*tinf);
mt2 = map_count_moment(mmap,t2,1:3);
m3t2 = mt2(3) - 3*mt2(2)*mt2(1) + 2*mt2(1)^3;

% per-class rates
ai = mmap_count_mean(mmap,1);

if strcmp(method,'exact_delta') == 1 || strcmp(method,'approx_delta') == 1
    % per-class variance difference
    dvt3 = zeros(m,1);
    for i = 1:m
        mmap2 = {mmap{1},mmap{2},mmap{2+i},mmap{2}-mmap{2+i}};
        Vt3 = mmap_count_var(mmap2,t3);
        dvt3(i) = Vt3(1)-Vt3(2);
    end
end

if strcmp(method,'exact_delta') == 1
    fit = m3pp2m_fit_count(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3);
elseif strcmp(method,'approx_delta') == 1
    fit = m3pp2m_fit_count_approx(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3);
elseif strcmp(method,'approx_cov') == 1
    vi = mmap_count_var(mmap,t3);
    s = 0.5 * (map_count_var(mmap,t3) - sum(vi));
    fit = m3pp22_fit_count_approx_cov(a, bt1, bt2, binf, m3t2, t1, t2, ai, s, t3);
elseif strcmp(method,'approx_ag') == 1
    gt3 = zeros(m,1);
    for i = 1:m
        mmap2 = {mmap{1},mmap{2},mmap{2+i},mmap{2}-mmap{2+i}};
        v2t3 = mmap_count_var(mmap2,t3);
        s2t3 = mmap_count_mcov(mmap2,t3);
        gt3(i) = v2t3(1) + s2t3(1,2);
    end
    fit = m3pp2m_fit_count_approx_ag(a, bt1, bt2, binf, m3t2, t1, t2, ai, gt3, t3);
else
    error('Invalid method ''%s\''', method);
end

end