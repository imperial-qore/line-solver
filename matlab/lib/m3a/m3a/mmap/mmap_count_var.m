function vk = mmap_count_var(mmap,t)
% Computes the variance of the counting process, at resolution t, for the
% given Marked MAP.
% Input:
% - mmap: the Marked MAP
% - t: the period considered for each sample of the counting process
% Output:
% - vk: the vector with the variance for each job class

n = size(mmap{1},1);
K = size(mmap,2)-2;

if map_issym(mmap)
    I = sym(eye(n));
    e = sym(ones(n,1));
    lk = sym(zeros(K,1));
    ck = sym(zeros(K,n));
    dk = sym(zeros(n,K));
    llk = sym(zeros(K,1));
    vk = sym(zeros(K,1));
else
    I = eye(n);
    e = ones(n,1);
    lk = zeros(K,1);
    ck = zeros(K,n);
    dk = zeros(n,K);
    llk = zeros(K,1);
    vk = zeros(K,1);
end

D = mmap{1} + mmap{2};
theta = map_prob(mmap);
tmp = (e * theta - D)^(-1);

for k=1:K
    % arrival rate
    lk(k) = theta * mmap{2+k} * e;
    ck(k,:) = theta * mmap{2+k} * tmp;
    dk(:,k) = tmp * mmap{2+k} * e;
    llk(k) = theta * mmap{2+k} * e;
end

for k=1:K
    vk(k) = (llk(k)-2*lk(k)^2 + 2*ck(k,:)*mmap{2+k}*e)*t - 2*ck(k,:)*(I-expm(D*t))*dk(:,k);
end

end