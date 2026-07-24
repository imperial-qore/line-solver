function MMAP = mamap2m_fit_gamma_fb_trace(T, A)
% Performs approximate fitting of a marked trace, yielding a second-order
% acyclic MMAP that matches the class probabilities, the forward and
% backward moments.
%
% Input
% - T: the inter-arrival times
% - A: the class marks of each job
% Output
% - MMAP: fitted MAMAP[m]

M1 = mean(T);
M2 = mean(T.^2);
M3 = mean(T.^3);
GAMMA = trace_gamma(T);

P = mtrace_pc(T,A);
F = mtrace_forward_moment(T,A,1);
B = mtrace_backward_moment(T,A,1);

MMAP = mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B);

end