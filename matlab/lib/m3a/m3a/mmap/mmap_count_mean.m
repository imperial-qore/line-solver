function [mk] = mmap_count_mean(mmap,t)
% Computes the mean of the counting process, at resolution t, for the
% given Marked MAP.
% Input:
% - mmap: the Marked MAP
% - t: the period considered for each sample of the counting process
% Output:
% - mk: the vector with the mean for each job class

if isempty(mmap)
    mk=[];
    return
end
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

% mean
mk = lk * t;

end