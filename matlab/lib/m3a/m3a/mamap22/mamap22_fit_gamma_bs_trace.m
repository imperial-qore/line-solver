function mmap = mamap22_fit_gamma_bs_trace(T, A)
% Performs approximate fitting of a marked trace with two classes, yielding
% a second-order acyclic MMAP[2] that fits the backward moments and the
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

P = mtrace_pc(T,A);
B = mtrace_backward_moment(T,A,1);
S = mtrace_sigma(T,A);

mmap = mamap22_fit_gamma_bs(M1, M2, M3, GAMMA, P, B, S);

end