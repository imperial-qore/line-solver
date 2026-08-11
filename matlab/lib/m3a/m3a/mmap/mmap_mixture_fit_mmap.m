function [MMAP,PHs] = mmap_mixture_fit_mmap(mmap)
% Fits a MMAP with m classes using a mixture of m^2 PH-distributions.
% Each PH distribution represents the probability distribution conditioned
% on the fact that the last arrival was of class i and the next arrival is
% of class j.
%
% INPUT
% - mmap: the MMAP to fit
% OUTPUT
% - MMAP: the fitted MMAP
% - PHs: the fitted PH-distributions for each transition

% compute two-step transition probabilities
P2 = mmap_sigma2(mmap);

% compute cross moments of order 1, 2 and 3
M1 = mmap_cross_moment(mmap,1);
M2 = mmap_cross_moment(mmap,2);
M3 = mmap_cross_moment(mmap,3);

% apply fitting algorithm
[MMAP,PHs] = mmap_mixture_fit(P2, M1, M2, M3);

end % end function