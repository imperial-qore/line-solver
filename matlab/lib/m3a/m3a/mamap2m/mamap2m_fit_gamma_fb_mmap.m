function MMAP = mamap2m_fit_gamma_fb_mmap(mmap)
% Performs approximate fitting of an MMAP[m], yielding a second-order
% acyclic MMAP[m] that matches the class probabilities, the forward and
% backward moments.
%
% Input
% - mmap: the MMAP[m] to fit
% Output
% - MMAP: fitted second-order MAMAP[m]

M1 = map_moment(mmap, 1);
M2 = map_moment(mmap, 2);
M3 = map_moment(mmap, 3);
GAMMA = map_gamma(mmap);

P = mmap_pc(mmap);
F = mmap_forward_moment(mmap,1);
B = mmap_backward_moment(mmap,1);

MMAP = mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B);

end