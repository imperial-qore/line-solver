function MMAP = mamap2m_fit_trace(T, C, fbsWeights)
% Fits a MAPH(2,m) or A MAMAP(2,m) that matches the characteristics of the
% trace. Three characteristics of the marked process are matched either
% exactly or approximately. Among these, the class probabilities are always
% matched exactly. The remaining two characteristics, by default, are the
% forward and backward moments, unless the underlying AMAP(2) is
% degenerate. Different pairs of characteristics can be chosen by
% specifying higher weights in the optional parameter fbsWeights.
%
% Input
% - T: the inter-arrival times
% - C: the class labels
% - fbsWeights: the weight assigned to forward moments, backward moments
%               and class transition probabilities (default: [1, 1, 1])
% Output
% - MMAP: fitted MAMAP[m]

if nargin == 2
    % by default equal weights
    % when weights are equal, F+B is preferred over F+S and B+S
    fbsWeights = [1 1 1];
end

M1 = mean(T);
M2 = mean(T.^2);
M3 = mean(T.^3);
GAMMA = trace_gamma(T);

P = mtrace_pc(T,C);
F = mtrace_forward_moment(T,C,1);
B = mtrace_backward_moment(T,C,1);
S = mtrace_sigma(T,C);

MMAP = mamap2m_fit(M1, M2, M3, GAMMA, P, F, B, S, fbsWeights);

end