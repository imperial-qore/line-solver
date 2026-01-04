function C=mtrace_iat2counts(T,A,scale,no_mex)
% Computes the per-class counting processes of T, i.e., the counts after
% "scale" units of time from an arrival.
% INPUT:
% - T: inter-arrival times
% - A: class labels
% - scale: time after an arrival
% - no_mex: set to 1 to avoid using the mex file (optional, default = 0)
% OUTPUT:
% - C: column k is the counting process for class k

if nargin < 4
    no_mex = 0;
end

if ~no_mex && exist('mtrace_iat2counts_native') == 3
   C = mtrace_iat2counts_native(T,A,scale);
else
    warning('Using slow Matlab function instead of MEX file.');
    n = length(T);
    CT = cumsum(T);
    K = unique(A); % classes
    C = zeros(n-1,length(K));
    for i=1:n-1
        if i >= 2
            % speedup loop by looking at the previous value
            cur = (i-1) + C(i-1);
        else
            cur = 1;
        end
        while CT(cur + 1) - CT(i) <= scale
            cur = cur + 1;
            % when the window first hits the end of the trace we return
            if cur == n 
                for j = 1:length(K)
                    C(i,j) = sum(A((i+1):cur) == K(j));
                end         
                C=C(1:i,:);
                return;
            end
        end
        for j = 1:length(K)
            C(i,j) = sum(A((i+1):cur) == K(j));
        end
    end
end

end
