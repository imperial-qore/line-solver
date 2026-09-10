function [fit,m3pps] = m3pp_superpos_fitc(av, btv, binfv, ...
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
    mmpp = mmpp2_fitc(av(i), ...
                           btv(i), btv(i), binfv(i), ...
                           m3tv(i), t, tinf);
    m3pps{i} = {mmpp{1},mmpp{2},mmpp{2}};
end

% perform superposition
fit = mmap_super(m3pps);

end