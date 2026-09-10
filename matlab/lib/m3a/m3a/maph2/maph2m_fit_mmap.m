function MAPH = maph2m_fit_mmap(mmap)
% Fits an MMAP[m] with a second-order MAPH[m] that matches the class
% probabilities (always fitted exactly) and the backward moments.
%
% Input
% - mmap: MMAP to fit (of arbitrary order)
% Output
% - MAPH: fitted second-order MAPH[m]

M1 = map_moment(mmap,1);
M2 = map_moment(mmap,2);
M3 = map_moment(mmap,3);

P = mmap_pc(mmap);
B = mmap_backward_moment(mmap,1);

MAPH = maph2m_fit(M1, M2, M3, P, B);

end