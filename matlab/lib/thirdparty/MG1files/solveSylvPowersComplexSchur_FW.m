function Y = solveSylvPowersComplexSchur_FW2(A, B, C)
% X = solveSylvPowersComplexSchur_FW() solves the equation 
% sum_{j=1}^N B_j Y A^{j-1} = C
% A: m x m
% B = [B_1 B_2 ... B_N], B_j: n x n
% C = n x m

m = size(A,1);
n = size(B,1);
N = size(B,2)/n;

[U,T] = schur(A,'complex');

E = C*U;

Tvec = zeros(m,N*m);
Tvec(:,1:m) = eye(m);
Tvec(:,m+1:2*m) = T;
for i = 2:N-1
    Tvec(:,i*m+1:(i+1)*m) = Tvec(:,(i-1)*m+1:i*m)*T;
end

%Z: matrix coeff of X_k
X = zeros(n*m,1);
k=1;
Bold = B;
B = reshape(B, n^2, N);
Tvec = reshape(Tvec, m^2, N).';
epsilon = 10e-14;
while k <= m
        if k > 1
            W = reshape(B*Tvec(:,(k-1)*m+1:(k-1)*m+k-1), n,n*(k-1))*X(1:n*(k-1));
        else
            W = zeros(n,1);
        end
        Z = reshape(B*Tvec(:,(k-1)*m+k), n, n);
        X(n*(k-1)+1:n*k) = Z\(E(:,k) - W);
        k = k+1;
end
Y = real(reshape(X,n,m)*U');