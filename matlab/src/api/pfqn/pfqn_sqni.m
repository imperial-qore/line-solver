%{
%{
 % @file pfqn_sqni.m
 % @brief Single Queue Network Iteration (SQNI) approximate solver.
%}
%}

%{
%{
 % @brief Single Queue Network Iteration (SQNI) approximate solver.
 % @fn pfqn_sqni(N, L, Z)
 % @param N Population vector.
 % @param L Service demand vector.
 % @param Z Think time vector.
 % @return Q Mean queue lengths.
 % @return U Utilization.
 % @return X System throughput.
%}
%}
function [Q,U,X]=pfqn_sqni(N,L,Z)
queueIdx = 1;
C = size(L,2);
Nt = sum(N);
X = zeros(1,C);
if sum(N)<=0
    Q = zeros(2,C);
    U = zeros(2,C);
    X = zeros(1,C);
elseif sum(N)==1
    for r=1:C
        X(r) = N(r)/(Z(r)+L(r));
        U(queueIdx,r) = X(r)*L(r);
        Q(queueIdx,r) = X(r)*L(r);
    end
else
    for r=1:C
        Nr = N(r);
        Lr = L(r);
        Zr = Z(r);
        Nvec_1r = N; Nvec_1r(r) = Nvec_1r(r) - 1;
        Br = N./(Z+L+L.*(sum(N)-1-sum(Z.*Nvec_1r./(Z+L+L*(sum(N)-2))))) .* Z;
        Br = Lr * sum(Br(setdiff(1:C,r)));
        if Lr==0
            X(r) = Nr / Zr;
        else
            X(r) = (Zr - (Br^2 - 2*Br*Lr*Nt - 2*Br*Zr + Lr^2*Nt^2 + 2*Lr*Nt*Zr - 4*Nr*Lr*Zr + Zr^2)^(1/2) - Br + Lr*Nt)/(2*Lr*Zr);
        end
        U(queueIdx,r) = X(r)*L(r);
        Q(queueIdx,r) = N(r)-X(r)*Z(r);
    end
end
for r=1:C
    if Z(r)==0
        X(r) = N(r)/(L(r)*(1+sum(Q,2,'omitnan')));
        U(queueIdx,r) = X(r)*L(r);
        Q(queueIdx,r) = N(r)-X(r)*Z(r);
    end
end
end