function MMAP = mmap_exponential(lambda, n)
% fits a order-n MMAP with given arrival rates lambda
if nargin < 2
    n = 1;
end
K = length(lambda);
MMAP = cell(1,2+K);
MMAP{1} = zeros(n);
MMAP{2} = zeros(n);
for k=1:K
    MMAP{2+k} = flip(eye(n))*lambda(k);
    MMAP{2} = MMAP{2} + MMAP{2+k};
end
MMAP = mmap_normalize(MMAP);
end