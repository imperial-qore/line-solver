function MAPH = maph2m_fit_trace(T,C)
% Fits a marked trace with a second-order MAPH[m] that matches the class
% probabilities (always fitted exactly) and the backward moments.
%
% Input
% - T: inter-arrival times
% - C: class labels
% - MAPH: fitted second-order MAPH[m]

M1 = mean(T);
M2 = mean(T.^2);
M3 = mean(T.^3);

P = mtrace_pc(T,C);
B = mtrace_backward_moment(T,C,1);

MAPH = maph2m_fit(M1, M2, M3, P, B);

end