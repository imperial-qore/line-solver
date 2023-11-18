function [fit,m3pps] = m3pp_superpos_fit_count(av, btv, binfv, ...
                                               m3tv, t, tinf)
% Fits k second-order M3PP[m_j] and superposes them into a
% M3PP[m] of order k+1, with m = \sum_j=1^k m_j.
%
% INPUT
%  av:      vector of length m with the per-process rates
%  btv:     vector of length k with the per-process IDC(t)
%  binfv:   vector of length k with the per-process IDC(inf)
%  m3tv:    third moment of counts
%  t:       finite time scale
%  tinf:    near-infinite time scale
% 
% OUTPUT
%  fit:     result of the superposition, M3PP[m] of order k+1
%  m3pps:   cell-array with the fitted and superposed second-order
%           M3PP[m_j]

% number of classes
m = length(av);

% total rate
a = sum(av);

% fit m3pp[2] processes
m3pps = cell(m,1);
for i = 1:m
    mmpp = mmpp2_fit_count(av(i), ...
                           btv(i), btv(i), binfv(i), ...
                           m3tv(i), t, tinf);
    m3pps{i} = {mmpp{1},mmpp{2},mmpp{2}};
end

% perform superposition
fit = mmap_superpos(m3pps);
% compare
fa = map_count_mean(fit,1);
fprintf('Rate: input = %f, %f\n', a, fa);
fav = mmap_count_mean(fit,1);
for i = 1:m
    fprintf('Class %d, rate: input = %f, output = %f\n', i, av(i), fav(i));
end
fbtv = mmap_count_var(fit,t) ./ (av * t);
for i = 1:m
    fprintf('Class %d, IDC(%.2f): input = %f, output = %f\n', i, t, btv(i), fbtv(i));
end
fbinfv = mmap_count_var(fit,tinf) ./ (av * tinf);
fbinfv_limit = mmap_count_var(fit, 1e6) ./ (av * 1e6);
for i = 1:m
    fprintf('Class %d, IDC(inf = %.2f): input = %f, output = %f, output_limit = %f\n', i, tinf, binfv(i), fbinfv(i), fbinfv_limit(i));
end

end