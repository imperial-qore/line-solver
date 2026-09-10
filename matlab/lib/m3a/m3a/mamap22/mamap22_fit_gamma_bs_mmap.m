function MMAP = mamap22_fit_gamma_bs_mmap(mmap)
% Performs approximate fitting of an MMAP[2], yielding a second-order
% acyclic MMAP[2] fitting the backward moments and the class transition
% probabilities.
%
% Input
% - mmap: the MMAP[2] to fit (arbitrary order)
% Output
% - MMAP: fitted second-order MAMAP[2]

M1 = map_moment(mmap,1);
M2 = map_moment(mmap,2);
M3 = map_moment(mmap,3);
GAMMA = map_gamma(mmap);

P = mmap_pc(mmap);
B = mmap_backward_moment(mmap,1);
S = mmap_sigma(mmap);

MMAP = mamap22_fit_gamma_bs(M1, M2, M3, GAMMA, P, B, S);

end