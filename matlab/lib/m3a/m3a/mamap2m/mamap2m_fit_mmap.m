function MMAP = mamap2m_fit_mmap(mmap, fbsWeights)
% Fits a MAPH(2,m) or A MAMAP(2,m) that matches the characteristics of the
% mmap. Three characteristics of the marked process are matched either
% exactly or approximately. Among these, the class probabilities are always
% matched exactly. The remaining two characteristics, by default, are the
% forward and backward moments, unless the underlying AMAP(2) is
% degenerate. Different pairs of characteristics can be chosen by
% specifying higher weights in the optional parameter fbsWeights.
%
% Input
% - mmap: the MMAP to fit
% - fbsWeights: the weight assigned to forward moments, backward moments
%               and class transition probabilities (default: [1, 1, 1])
% Output
% - MMAP: fitted MAMAP[m]

if nargin == 1
    % by default equal weights
    % when weights are equal, F+B is preferred over F+S and B+S
    fbsWeights = [1 1 1];
end

M1 = map_moment(mmap,1);
M2 = map_moment(mmap,2);
M3 = map_moment(mmap,3);
GAMMA = map_gamma(mmap);

P = mmap_pc(mmap);
F = mmap_forward_moment(mmap,1);
B = mmap_backward_moment(mmap,1);
S = mmap_sigma(mmap);

MMAP = mamap2m_fit(M1, M2, M3, GAMMA, P, F, B, S, fbsWeights);

end