function [FIT] = m3pp22_fit_count_approx_cov(a, bt1, bt2, binf, m3t2, ...
                                             t1, t2, ...
                                             ai, st3, t3)
% Fits a second-order Marked MMPP.
% INPUT:
% - a: arrival rate
% - bt1: IDC at scale t1
% - bt2: IDC at scale t2
% - binf: IDC for t->inf
% - m3t2: third central moment
% - t1: first time scale
% - t2: second time scale
% - ai: rates of the two classes
% - st3: count covariance between the two classes at scale t3
% - t3: third time scale

if abs(a - sum(ai)) > 1e-8
    error('Inconsistent per-class arrival rates.');
end

% fit underlying MMPP(2)
FIT = mmpp2_fit_count_approx(a, bt1, bt2, binf, m3t2, t1, t2);

% fit M3PP(2,2)
FIT = m3pp22_fit_count_approx_cov_multiclass(FIT, ai, st3, t3);


end