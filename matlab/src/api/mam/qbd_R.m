function R=qbd_R(B,L,F,iter_max)
% Successive substitutions method
if nargin<4
    iter_max = 100000;
end
Fil = F * inv(L);
BiL = B * inv(L);
R = -Fil;
Rprime = -Fil - R^2*BiL;
for iter=1:iter_max
    R = Rprime;
    Rprime = -Fil -R^2*BiL;
    if norm(R-Rprime,1)<=1e-12
        break
    end
end
R = Rprime;
end