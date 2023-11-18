function [fit] = m3pp2m_fit_count_trace(T, A, method, t1, tinf)
% Fits the theoretical characteristics of a MMAP(n,m) with a M3PP(2,m).
% INPUT:
% - mmap: the MMAP(n,m) to fit with a M3PP(2,m)
% - method: 'exact_delta', 'approx_delta', 'approx_cov' or 'approx_ag'

labels = unique(A);
m = length(labels);

if nargin < 3
    method = 'approx_delta';
end

if strcmp(method,'approx_cov') && m > 2
   error('Approximate covariance fitting only supported for two classes.');
end

TC = cumsum(T);

if nargin <= 3
    t1 = 10*mean(T);
    tinf = max(10*t1, (TC(end)-TC(1)) / 100);
end
t2 = t1+mean(T);
t3 = tinf; % this controls the approximation - GC changed from t1 to tinf

fprintf('Computing counting process at resolution %f\n', t1);
mNt1 = mtrace_iat2counts(T, A, t1);
mNt2 = mNt1;
fprintf('Computing counting process at resolution %f\n', tinf);
mNtinf = mtrace_iat2counts(T, A, tinf);
mNt3 = mNt1;

Nt1 = sum(mNt1, 2);
Nt2 = sum(mNt2, 2);
Ntinf = sum(mNtinf, 2);

% total rate
a = 1/mean(T);

% per-class rates
ai = zeros(m,1);
for i = 1:m
    ai(i) = a * sum(A == labels(i))/length(A);
end

fprintf('Rate: %f\n', a);

% joint-process charactersitics
bt1 = var(Nt1)/(a*t1);
bt2 = bt1;
binf = var(Ntinf)/(a*tinf);
mt2 = [mean(Nt2), mean(Nt2.^2), mean(Nt2.^3)];
m3t2 = mt2(3) - 3*mt2(2)*mt2(1) + 2*mt2(1)^3;

if strcmp(method,'exact_delta') == 1 || strcmp(method,'approx_delta') == 1
    % per-class variance differential
    dvt3 = zeros(m,1);
    for i = 1:m
        N2 = [mNt3(:,i), sum(mNt3,2)-mNt3(:,i)];
        Vt3 = [var(N2(:,1)), var(N2(:,2))];
        dvt3(i) = Vt3(1)-Vt3(2);
    end
end

if strcmp(method,'exact_delta') == 1
    fit = m3pp2m_fit_count(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3);
elseif strcmp(method,'approx_delta') == 1
    fit = m3pp2m_fit_count_approx(a, bt1, bt2, binf, m3t2, t1, t2, ai, dvt3, t3);
elseif strcmp(method,'approx_cov') == 1
    V = var(sum(mNt3,2));
    vi = [var(mNt3(:,1)), var(mNt3(:,2))];
    s = 0.5 * (V - sum(vi));
    fit = m3pp22_fit_count_approx_cov(a, bt1, bt2, binf, m3t2, t1, t2, ai, s, t3);
elseif strcmp(method,'approx_ag') == 1
    vt3 = var(mNt3)';
    st3 = zeros(m,1);
    for i = 1:m
        covt3 = cov(mNt3(:,i), sum(mNt3,2)-mNt3(:,i));
        st3(i) = covt3(1,2);
    end
    gt3 = vt3 + st3;
    fit = m3pp2m_fit_count_approx_ag(a, bt1, bt2, binf, m3t2, t1, t2, ai, gt3, t3);
else
    error('Invalid method ''%s\''', method);
end

end