function [lk] = mmap_count_lambda(mmap)
% Computes the arrival rate of the counting process, for the
% given Marked MAP.
% Input:
% - mmap: the Marked MAP
% Output:
% - lk: the vector with the rate for each job class

n = size(mmap{1},1);
K = size(mmap,2)-2;

if (map_issym(mmap))
    e = sym(ones(n,1));
    lk = sym(zeros(K,1));
else
    e = ones(n,1);
    lk = zeros(K,1);
end

theta = map_prob(mmap);

for k=1:K
    lk(k) = theta * mmap{2+k} * e;
end
lk=lk';
end