function Y = solveSylvPowersDirectSum(A, B, C)
% X = solveSylvPowersDirectSum() solves the equation 
% sum_{j=1}^N B_j Y A^{j-1} = C, directly
% A: m x m
% B = [B_1 B_2 ... B_N], B_j: n x n
% C = n x m


m = size(A,1);
n = size(B,1);
N = size(B,2)/n;


Avec = zeros(m,N*m);
Avec(:,1:m) = eye(m);
Avec(:,m+1:2*m) = A';
for i = 2:N-1
    Avec(:,i*m+1:(i+1)*m) = Avec(:,(i-1)*m+1:i*m)*Avec(:,m+1:2*m);
end

%Z: matrix coeff of vec(Y)
Z = zeros(m*n);
for j = 1:N
    Z = Z + kron(Avec(:,(j-1)*m+1:j*m), B(:,(j-1)*n+1:j*n));
end

Y = Z\(reshape(C,m*n,1));

Y = reshape(Y, n, m);