function mmap = mamap22_fit_gamma_fs_trace(T, C)
% Performs approximate fitting of a marked trace with two classes, yielding
% a second-order acyclic MMAP[2] that fits the forward moments and the
% class transition probabilities.
%
% Input
% - T: the inter-arrival times
% - C: the class marks of each job
% Output
% - mmap: fitted MAMAP[m]

M1 = mean(T);
M2 = mean(T.^2);
M3 = mean(T.^3);
GAMMA = trace_gamma(T);

P = mtrace_pc(T,C);
F = mtrace_forward_moment(T,C,1);
S = mtrace_sigma(T,C);

mmap = mamap22_fit_gamma_fs(M1, M2, M3, GAMMA, P, F, S);

end