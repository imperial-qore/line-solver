function summary = mtrace_summary(T,C)

summary = struct();
summary.M = [mean(T), mean(T.^2), mean(T.^3), mean(T.^4), mean(T.^5)];
summary.ACF = trace_acf(T, 1:100);
summary.F1 = mtrace_forward_moment(T, C, 1);
summary.F2 = mtrace_forward_moment(T, C, 2);
summary.B1 = mtrace_backward_moment(T, C, 1);
summary.B2 = mtrace_backward_moment(T, C, 2);
summary.C1 = mtrace_cross_moment(T, C, 1);
summary.C2 = mtrace_cross_moment(T, C, 2);
summary.Pc = mtrace_pc(T, C);
summary.Pab = mtrace_sigma(T, C);