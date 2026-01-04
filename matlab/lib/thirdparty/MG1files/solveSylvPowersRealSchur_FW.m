function Y = solveSylvPowersRealSchur_FW(A, B, C)
% X = solveSylvPowersRealSchur_FW() solves the equation 
% sum_{j=1}^N B_j Y A^{j-1} = C
% A: m x m
% B = [B_1 B_2 ... B_N], B_j: n x n
% C = n x m

m = size(A,1);
n = size(B,1);
N = size(B,2)/n;

[U,T] = schur(A,'real');

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
Tvec = reshape(Tvec, m^2, N)';
epsilon = 10e-14;
while k <=m-1
    if abs(T(k+1,k)) < epsilon
        if k > 1
            W = reshape(B*Tvec(:,(k-1)*m+1:(k-1)*m+k-1), n,n*(k-1))*X(1:n*(k-1));
        else
            W = zeros(n,1);
        end
        Z = reshape(B*Tvec(:,(k-1)*m+k), n, n);
        
        X(n*(k-1)+1:n*k) = Z\(E(:,k) - W);
        k = k+1;
    else
        if k > 1
            W = [  reshape(B*Tvec(:,(k-1)*m+1:(k-1)*m+k-1), n, n*(k-1))*X(1:n*(k-1));
                   reshape(B*Tvec(:,k*m+1:k*m+k-1),         n, n*(k-1))*X(1:n*(k-1))];
        else
            W = zeros(2*n,1);
        end
        Z = [   reshape(B*Tvec(:,(k-1)*m+k), n, n)     reshape(B*Tvec(:,(k-1)*m+k+1), n, n);
                reshape(B*Tvec(:,k*m+k), n, n)   reshape(B*Tvec(:,k*m+k+1), n, n)];
        
        X(n*(k-1)+1:n*(k+1)) = Z\([E(:,k);E(:,k+1)] - W);
        k = k+2;
    end
end
      
if k == m
    W = reshape(B*Tvec(:,(k-1)*m+1:(k-1)*m+k-1), n, n*(k-1))*X(1:n*(k-1));
    Z = reshape(B*Tvec(:,(k-1)*m+k), n, n);
    X(n*(k-1)+1:n*k) = Z\(E(:,k) - W);
end

Y = reshape(X,n,m)*U';
