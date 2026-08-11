function P = mxpow(A, k)
%MXPOW  Matrix power returning the identity when k==0.
%   P = MXPOW(A, k) returns A^k, or eye(size(A,1)) when k==0. Used by the
%   transient-QBD solver where zero-length level ranges must yield identity.
if k == 0
    P = eye(size(A, 1));
else
    P = A^k;
end
end
