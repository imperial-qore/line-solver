function [FIT] = m3pp2m_fit_count_approx_ag(a, bt1, bt2, binf, m3t2, ...
                                            t1, t2, ...
                                            ai, gt3, t3)
% Fits a second-order Marked MMPP.
% a: arrival rate
% bt1: IDC at scale t1
% bt2: IDC at scale t2
% binf: IDC for t->inf
% m3t2: third central moment
% t1: first time scale
% t2: second time scale
% ai: i-th element is the rate of class i
% gt3: i-th element is the sum of the variance of class i and its
%      covariance with all the other classes combined.
% t3: third time scale

if abs(a - sum(ai)) > 1e-8
    error('Inconsistent per-class arrival rates.');
end

% fit underlying MMPP(2)
FIT = mmpp2_fit_count_approx(a, bt1, bt2, binf, m3t2, t1, t2);

% fit M3PP(2,m)
FIT = m3pp2m_fit_count_approx_ag_multiclass(FIT, ai, gt3, t3);

end