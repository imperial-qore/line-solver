function MMAP = mamap22_fit_gamma_fs_mmap(mmap)
% Performs approximate fitting of an MMAP[2], yielding a second-order
% acyclic MMAP[2] fitting the forward moments and the class transition
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
F = mmap_forward_moment(mmap,1);
S = mmap_sigma(mmap);

MMAP = mamap22_fit_gamma_fs(M1, M2, M3, GAMMA, P, F, S);

end