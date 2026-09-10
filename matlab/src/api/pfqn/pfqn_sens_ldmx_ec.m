%{
%{
 % @file pfqn_sens_ldmx_ec.m
 % @brief Effective capacity terms of the mixed load-dependent MVA together
 %        with their exact derivatives with respect to the open-class load.
%}
%}

%{
%{
 % @brief Computes the effective capacity terms EC, E and Eprime of the mixed
 %        load-dependent MVA of Bruell-Balbo-Afshari, exactly as pfqn_ldmx_ec
 %        does, and additionally their exact analytic derivatives with respect
 %        to the open-class load Lo(i) of each station.
 %
 %        Lo(i) = sum_r lambda(r)*D(i,r) is the only channel through which a
 %        service demand enters E, Eprime and EC: the load-dependent rates mu
 %        are independent of the demands. Station i's terms depend on Lo(i)
 %        alone, so a single derivative per station is enough, and the chain
 %        rule then yields the derivative with respect to any demand-scaling
 %        parameter. This is the factorization behind equations (19), (21) and
 %        (24)-(31) of the reference, which are reproduced here term by term.
 %
 %        Reference: I. F. Akyildiz and J. C. Strelen, "Moment Analysis for
 %        Load-Dependent Mixed Product Form Queueing Networks", IEEE Trans.
 %        Communications 39(6):828-832, 1991.
 %
 % @fn pfqn_sens_ldmx_ec(lambda, D, mu)
 % @param lambda Arrival rate vector (1 x R). Zero for closed classes.
 % @param D Service demand matrix (M x R).
 % @param mu Load-dependent rate matrix (M x Nt), limited load dependence.
 % @return EC Effective capacity matrix (M x Nt).
 % @return E E-function values (M x 1+Nt).
 % @return Eprime E-prime function values (M x 1+Nt).
 % @return Lo Open class load vector (M x 1).
 % @return dEC dEC(i,n)/dLo(i), same shape as EC.
 % @return dE dE(i,1+n)/dLo(i), same shape as E.
 % @return dEprime dEprime(i,1+n)/dLo(i), same shape as Eprime.
%}
%}
function [EC,E,Eprime,Lo,dEC,dE,dEprime] = pfqn_sens_ldmx_ec(lambda,D,mu)
% [EC,E,EPRIME,LO,DEC,DE,DEPRIME] = PFQN_SENS_LDMX_EC(LAMBDA,D,MU)
[M,~] = size(mu);
Lo = zeros(M,1);
for ist=1:M
    Lo(ist) = lambda*D(ist,:)';
end

b = zeros(M,1); % limited load dependence level
for ist=1:M
    b(ist) = find(mu(ist,:)==mu(ist,end), 1 );
end
Nt = size(mu,2);
mu(:,end+1:end+1+max(b)) = repmat(mu(:,end),1,1+max(b));
C = 1./mu;

EC = zeros(M,Nt);
E = zeros(M,1+Nt);
Eprime = zeros(M,1+Nt);
dEC = zeros(M,Nt);
dE = zeros(M,1+Nt);
dEprime = zeros(M,1+Nt);

for ist=1:M
    bi = b(ist);
    Cb = C(ist,bi);
    Lo_i = Lo(ist);
    den = 1 - Lo_i*Cb;   % geometric tail factor of the limited load dependence

    E1 = zeros(1,1+Nt);   dE1 = zeros(1,1+Nt);
    E2 = zeros(1,1+Nt);   dE2 = zeros(1,1+Nt);
    E3 = zeros(1,1+Nt);   dE3 = zeros(1,1+Nt);
    E2prime = zeros(1,1+Nt); dE2prime = zeros(1,1+Nt);
    F2 = zeros(1+Nt,1+max(0,bi-2));  dF2 = zeros(1+Nt,1+max(0,bi-2));
    F3 = zeros(1+Nt,1+max(0,bi-2));  dF3 = zeros(1+Nt,1+max(0,bi-2));
    F2prime = zeros(1+Nt,1+max(0,bi-2)); dF2prime = zeros(1+Nt,1+max(0,bi-2));

    for n=0:Nt
        if n >= bi
            % E(n) = 1/den^(n+1)  =>  dE/dLo = (n+1)*Cb/den^(n+2)
            E(ist,1+n) = 1 / den^(n+1);
            dE(ist,1+n) = (n+1)*Cb / den^(n+2);
            Eprime(ist,1+n) = Cb*E(ist,1+n);
            dEprime(ist,1+n) = Cb*dE(ist,1+n);
        else % n <= bi-1
            %% E1 and its derivative, eq. (25)-(26)
            if n==0
                E1(1+n) = 1 / den;
                dE1(1+n) = Cb / den^2;
                for j=1:(bi-1)
                    E1(1+n) = E1(1+n) * C(ist,j) / Cb;
                    dE1(1+n) = dE1(1+n) * C(ist,j) / Cb;
                end
            else % n>0
                fac = Cb / C(ist,n);
                E1(1+n) = (1/den) * fac * E1(1+(n-1));
                dE1(1+n) = (Cb/den^2) * fac * E1(1+(n-1)) + (1/den) * fac * dE1(1+(n-1));
            end

            %% F2 and its derivative, eq. (27)-(28)
            for n0 = 0:(bi-2)
                if n0 == 0
                    F2(1+n,1+n0) = 1;
                    dF2(1+n,1+n0) = 0;
                else
                    coef = (n+n0)/n0 * C(ist,n+n0);
                    F2(1+n,1+n0) = coef * Lo_i * F2(1+n,1+(n0-1));
                    dF2(1+n,1+n0) = coef * (F2(1+n,1+(n0-1)) + Lo_i*dF2(1+n,1+(n0-1)));
                end
            end
            E2(1+n) = sum(F2(1+n,1+(0:bi-2)));
            dE2(1+n) = sum(dF2(1+n,1+(0:bi-2)));

            %% F3 and its derivative, eq. (29)-(30)
            for n0 = 0:(bi-2)
                if n == 0 && n0 == 0
                    F3(1+n,1+n0) = 1;
                    for j=1:(bi-1)
                        F3(1+n,1+n0) = F3(1+n,1+n0) * C(ist,j) / Cb;
                    end
                    dF3(1+n,1+n0) = 0;
                elseif n > 0 && n0 == 0
                    fac = Cb / C(ist,n);
                    F3(1+n,1+n0) = fac * F3(1+(n-1),1+0);
                    dF3(1+n,1+n0) = fac * dF3(1+(n-1),1+0);
                else
                    coef = (n+n0)/n0 * Cb;
                    F3(1+n,1+n0) = coef * Lo_i * F3(1+n,1+(n0-1));
                    dF3(1+n,1+n0) = coef * (F3(1+n,1+(n0-1)) + Lo_i*dF3(1+n,1+(n0-1)));
                end
            end
            E3(1+n) = sum(F3(1+n,1+(0:bi-2)));
            dE3(1+n) = sum(dF3(1+n,1+(0:bi-2)));

            %% F2prime and its derivative
            for n0 = 0:(bi-2)
                if n0 == 0
                    F2prime(1+n,1+n0) = C(ist,n+1);
                    dF2prime(1+n,1+n0) = 0;
                else
                    coef = (n+n0)/n0 * C(ist,n+n0+1);
                    F2prime(1+n,1+n0) = coef * Lo_i * F2prime(1+n,1+(n0-1));
                    dF2prime(1+n,1+n0) = coef * (F2prime(1+n,1+(n0-1)) + Lo_i*dF2prime(1+n,1+(n0-1)));
                end
            end
            E2prime(1+n) = sum(F2prime(1+n,1+(0:(bi-2))));
            dE2prime(1+n) = sum(dF2prime(1+n,1+(0:(bi-2))));

            % E = E1 + E2 - E3, eq. (23)-(24)
            E(ist,1+n) = E1(1+n) + E2(1+n) - E3(1+n);
            dE(ist,1+n) = dE1(1+n) + dE2(1+n) - dE3(1+n);
            if n<bi-1
                Eprime(ist,1+n) = Cb * E1(1+n) + E2prime(1+n) - Cb * E3(1+n);
                dEprime(ist,1+n) = Cb * dE1(1+n) + dE2prime(1+n) - Cb * dE3(1+n);
            else %n>=bi-1
                Eprime(ist,1+n) = Cb * E(ist,1+n);
                dEprime(ist,1+n) = Cb * dE(ist,1+n);
            end
        end
    end

    % EC(n) = C(n)*E(n)/E(n-1), eq. (19); quotient rule for the derivative
    for n=1:Nt
        EC(ist,n) = C(ist,n) * E(ist,1+n) / E(ist,1+(n-1));
        dEC(ist,n) = C(ist,n) * (dE(ist,1+n)*E(ist,1+(n-1)) - E(ist,1+n)*dE(ist,1+(n-1))) ...
                     / E(ist,1+(n-1))^2;
    end
end
end
