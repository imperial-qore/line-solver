function [fit] = maph2m_fit_count_theoretical(mmap, method)
% Fits the theoretical characteristics of a MMAP(n,m) with a M3PP(2,m).
% INPUT:
% - mmap: the MMAP(n,m) to fit with a M3PP(2,m)
% - method: either 'exact' or 'approx'

if nargin == 1
    method = 'exact';
end

m = size(mmap,2)-2;

t1 = 1;
t2 = 10;
tinf = 1e4;
t3 = 10;

% joint-process charactersitics
a = map_count_mean(mmap,t1)/t1;
bt1 = map_count_var(mmap,t1)/(a*t1);
binf = map_count_var(mmap,tinf)/(a*tinf);

% per-class rates
ai = mmap_count_mean(mmap,1);

% per-class variance differential
dvt3 = zeros(m,1);
for i = 1:m
    mmap2 = {mmap{1},mmap{2},mmap{2+i},mmap{2}-mmap{2+i}};
    Vt3 = mmap_count_var(mmap2,t3);
    dvt3(i) = Vt3(1)-Vt3(2);
end

if strcmp(method,'exact') == 1
    fit = maph2m_fit_count(a, bt1, binf, t1, ai, dvt3, t3);
elseif strcmp(method,'approx') == 1
    fit = maph2m_fit_count_approx(a, bt1, binf, t1, ai, dvt3, t3);
else
    error('Invalid method ''%s\''', method);
end

end