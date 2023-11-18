function [ C ] = kronsum( A, B )

% return the kronecker sum of two matrices
C = kron(A,eye(size(B,2))) + kron(eye(size(A,2)),B);

end

