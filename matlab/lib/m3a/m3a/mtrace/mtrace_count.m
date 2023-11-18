function N = mtrace_count(T,A,t)
% Computes the count process sample, at resolution, for the marked trace
% (T,A).

time = sum(T(2:end));
periods = ceil(time / t);
C = unique(A);

N = zeros(periods,length(C));
Tcum = cumsum(T);

for i = 1:periods
    tstart = T(1) + (i-1) * t;
    tend = T(1) + i * t;
    for c = 1:length(C)
        N(i,c) = sum( Tcum > tstart & Tcum <= tend & A == c );
    end
end
