function [fit,m3pps] = m3pp_superpos_fit_count_theoretical(MMAP, t, tinf)
% Superposes k M3PP to fit the characteristics of a MMAP[k].
% INPUT
% - MMAP: process to fit
% - t: finite time scale
% - tinf: near-infinite time scale

% number of classes
m = size(MMAP,2)-2;

% per-class rates
av = mmap_count_mean(MMAP,1);

% compute per-class IDC(t), IDC(inf)
btv = mmap_count_idc(MMAP,t);
binfv = mmap_count_idc(MMAP,tinf);

% compute per-class third central moment
mtv = mmap_count_moment(MMAP,t,1:3);
m3tv = zeros(m,1);
for i = 1:m
    m3tv(i) = mtv(3,i) - 3*mtv(2,i)*mtv(1,i) + 2*mtv(1,i)^3;
end

% fit superposition
[fit, m3pps] = m3pp_superpos_fit_count(av, btv, binfv, m3tv, t, tinf);

end