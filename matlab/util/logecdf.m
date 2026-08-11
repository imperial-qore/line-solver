function logecdf(S)
% LOGECDF(S)
% Plot an empirical cumulative distribution in log scale for dataset S
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[F,X]=ecdf(S);
semilogx(X,F);
xlabel('x');
ylabel('F(x)');
end