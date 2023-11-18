function R=qbd_R_logred(B,L,F,iter_max)
% Logarithmic reduction method
if nargin<4
    iter_max = 100000;
end
iLF = -inv(L)*F;
iLB = -inv(L)*B;
T = iLF;
S = iLB;
for iter=1:iter_max
    D = iLF*iLB + iLB*iLF;
    iLF = inv(eye(r)-D) *iLF*iLF;
    iLB = inv(eye(r)-D) *iLB*iLB;
    S = S + T*iLB;
    T = T*iLF;
    if norm(ones(r,1)-S*ones(r,1) ,1) <= 1e-12
        break
    end
end
U = L + F*S;
R = -F * inv(U);
end