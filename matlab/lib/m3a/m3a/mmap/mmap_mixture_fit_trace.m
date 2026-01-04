function [MMAP,PHs] = mmap_mixture_fit_trace(T, A)
% Fits a MMAP with m classes using a mixture of m^2 PH-distributions.
% Each PH distribution represents the probability distribution conditioned
% on the fact that the last arrival was of class i and the next arrival is
% of class j.
%
% INPUT
% - T: inter-arrival times
% - A: class labels
% OUTPUT
% - MMAP: the fitted MMAP
% - PHs: the fitted PH-distributions for each transition

% compute two-step transition probabilities
P2 = mtrace_sigma2(T,A);

% compute cross moments of order 1, 2 and 3
M1 = mtrace_cross_moment(T,A,1);
M2 = mtrace_cross_moment(T,A,2);
M3 = mtrace_cross_moment(T,A,3);

% apply fitting algorithm
[MMAP,PHs] = mmap_mixture_fit(P2, M1, M2, M3);

end % end function