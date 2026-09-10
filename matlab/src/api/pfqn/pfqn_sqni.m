%{
%{
 % @file pfqn_sqni.m
 % @brief Square-root Non-iterative (SQNI) approximate solver.
%}
%}

%{
%{
 % @brief Square-root Non-iterative (SQNI) approximate solver.
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
Q = zeros(1,C);
U = zeros(1,C);
if sum(N)<=0
    return
elseif sum(N)==1
    for r=1:C
        X(r) = N(r)/(Z(r)+L(r));
        U(queueIdx,r) = X(r)*L(r);
        Q(queueIdx,r) = X(r)*L(r);
    end
else
    % Pre-set Q for Z=0 chains (e.g., self-looping classes) to avoid NaN
    for r=1:C
        if Z(r)==0
            Q(queueIdx,r) = N(r);
        end
    end
    for r=1:C
        if Z(r)==0
            continue % handled after the main loop
        end
        Nr = N(r);
        Lr = L(r);
        Zr = Z(r);
        Nvec_1r = N;
        Nvec_1r(r) = Nvec_1r(r) - 1;
        Br = N./(Z+L+L.*(sum(N)-1-sum(Z.*Nvec_1r./(Z+L+L*(sum(N)-2))))) .* Z;
        Brsum = 0;
        for c = 1:C
            if c ~= r
                Brsum = Brsum + Br(c);
            end
        end
        Br = Lr * Brsum;
        if Lr==0
            X(r) = Nr / Zr;
        else
            discriminant = Br^2 - 2*Br*Lr*Nt - 2*Br*Zr + Lr^2*Nt^2 + 2*Lr*Nt*Zr - 4*Nr*Lr*Zr + Zr^2;
            if discriminant < 0
                discriminant = 0;
            end
            X(r) = (Zr - sqrt(discriminant) - Br + Lr*Nt)/(2*Lr*Zr);
        end
        U(queueIdx,r) = X(r)*L(r);
        Q(queueIdx,r) = N(r)-X(r)*Z(r);
    end
end
for r=1:C
    if Z(r)==0
        X(r) = N(r)/(L(r)*(1+sum(Q,2)));
        U(queueIdx,r) = X(r)*L(r);
        Q(queueIdx,r) = N(r)-X(r)*Z(r);
    end
end
end
