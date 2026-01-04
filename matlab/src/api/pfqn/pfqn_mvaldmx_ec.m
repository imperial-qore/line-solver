%{
%{
 % @file pfqn_mvaldmx_ec.m
 % @brief Compute effective capacity terms for MVALDMX solver.
%}
%}

%{
%{
 % @brief Compute effective capacity terms for MVALDMX solver.
 % @fn pfqn_mvaldmx_ec(lambda, D, mu)
 % @param lambda Arrival rate vector.
 % @param D Service demand matrix.
 % @param mu Load-dependent rate matrix.
 % @return EC Effective capacity matrix.
 % @return E E-function values.
 % @return Eprime E-prime function values.
 % @return Lo Open class load vector.
%}
%}
function [EC,E,Eprime,Lo] = pfqn_mvaldmx_ec(lambda,D,mu)
% [EC,E,EPRIME,LO] = PFQN_MVALDMX_EC(LAMBDA,D,MU)
% Compute the effective capacity terms in MVALDMX
% Think times are not handled since this assumes limited load-dependence
[M,R] = size(mu);
%Nt = sum(N(isfinite(N)));
Lo = zeros(M,1);
for ist=1:M
    Lo(ist) = lambda*D(ist,:)';
end

b = zeros(M,1); % limited load dependence level
for ist=1:M
    b(ist) = find(mu(ist,:)==mu(ist,end), 1 );
end
Nt = size(mu,2); % compute extra elements if present
mu(:,end+1:end+1+max(b)) = repmat(mu(:,end),1,1+max(b));
C = 1./mu;

EC = zeros(M,Nt);
%Ever = zeros(M,1+Nt);
E = zeros(M,1+Nt);
Eprime = zeros(M,1+Nt);
for ist=1:M
    E1 = zeros(1+Nt);
    E2 = zeros(1+Nt);
    E3 = zeros(1+Nt);
    F2 = zeros(1+Nt,1+b(ist)-2);
    F3 = zeros(1+Nt,1+b(ist)-2);
    
    E2prime = zeros(1+Nt);
    F2prime = zeros(1+Nt,1+b(ist)-2);
    for n=0:Nt
        if n >= b(ist)
            E(ist,1+n) = 1 / (1-Lo(ist)*C(ist,b(ist)))^(n+1);
            %            Ever(i,1+n) = 1 / (1-Lo(i)*C(i,b(i)))^(n+1);
            Eprime(ist,1+n) = C(ist,b(ist))*E(ist,1+n);
        else % n <= b(i)-1
            %% compute E1
            if n==0
                E1(1+n) = 1 / (1-Lo(ist)*C(ist,b(ist)));
                for j=1:(b(ist)-1)
                    E1(1+n) = E1(1+n) * C(ist,j) / C(ist,b(ist));
                end
            else % n>0
                E1(1+n) = 1 / (1-Lo(ist)*C(ist,b(ist))) * C(ist,b(ist)) / C(ist,n) * E1(1+(n-1));
            end
            
            %% compute F2
            for n0 = 0:(b(ist)-2)
                if n0 == 0
                    F2(1+n,1+n0) = 1;
                else
                    F2(1+n,1+n0) = (n+n0)/n0 * Lo(ist) * C(ist,n+n0) * F2(1+n,1+(n0-1));
                end
            end
            
            %% compute E2
            E2(1+n) = sum(F2(1+n,1+(0:b(ist)-2)));
            
            %% compute F3
            for n0 = 0:(b(ist)-2)
                if n == 0 && n0 == 0
                    F3(1+n,1+n0) = 1;
                    for j=1:(b(ist)-1)
                        F3(1+n,1+n0) = F3(1+n,1+n0) * C(ist,j) / C(ist,b(ist));
                    end
                elseif n > 0 && n0 == 0
                    F3(1+n,1+n0) = C(ist,b(ist)) / C(ist,n) * F3(1+(n-1),1+0);
                else
                    F3(1+n,1+n0) = (n+n0)/n0 * Lo(ist) * C(ist,b(ist)) * F3(1+n,1+(n0-1));
                end
            end
            
            %% compute E3
            E3(1+n) = sum(F3(1+n,1+(0:b(ist)-2)));
            
            %% compute F2prime
            for n0 = 0:(b(ist)-2)
                if n0 == 0
                    F2prime(1+n,1+n0) = C(ist,n+1);
                else
                    F2prime(1+n,1+n0) = (n+n0)/n0 * Lo(ist) * C(ist,n+n0+1) * F2prime(1+n,1+(n0-1));
                end
            end
            
            %% compute E2prime
            E2prime(1+n) = sum(F2prime(1+n,1+(0:(b(ist)-2))));
            
            % finally, compute E, Eprime, and EC
            E(ist,1+n) = E1(1+n) + E2(1+n) - E3(1+n);
            if n<b(ist)-1
                Eprime(ist,1+n) = C(ist,b(ist)) * E1(1+n) + E2prime(1+n) - C(ist,b(ist)) * E3(1+n);
            else %n>=b(i)-1
                Eprime(ist,1+n) = C(ist,b(ist)) * E(ist,1+n);
            end
            %% verification of E
            %            Ever(i,1+n) = C(i,b(i))^(n+1-b(i)) / (1-Lo(i)*C(i,b(i)))^(n+1) * prod(C(i,(n+1):(b(i)-1)));
            %            for n0 = 0:(b(i)-2)
            %                Ever(i,1+n) = Ever(i,1+n) + nchoosek(n+n0,n0) * Lo(i)^n0 * (prod(C(i,(n+1):(n+n0)))-C(i,b(i))^(n0+n+1-b(i))*prod(C(i,(n+1):(b(i)-1))));
            %            end
        end
        %         %% verification2 of E
        %         Ever(i,1+n) = 0;
        %         for n0=0:1000
        %             if n+n0+1>b(i)
        %                 C(i,n+n0+1) = 1/mu(i,b(i));
        %             end
        %             Ever(i,1+n) = Ever(i,1+n)  +nchoosek(n+n0,n0) * Lo(i)^n0 * prod(C(i,(n+1):(n+n0)));
        %         end
        %         %% verification of Eprime
        %         Eprimever(i,1+n) = 0;
        %         for n0=0:1000
        %             if n+n0+1>b(i)
        %                 C(i,n+n0+1) = 1/mu(i,b(i));
        %             end
        %             Eprimever(i,1+n) = Eprimever(i,1+n)  +nchoosek(n+n0,n0) * Lo(i)^n0 * prod(C(i,(n+1):(n+n0+1)));
        %         end
    end
    %    Eprime = Eprimever;
    % EC not defined for n=0
    for n=1:Nt
        %if n>=b(i)
        %    EC(i,n) = C(i,b(i)) / (1-Lo(i)*C(i,b(i)));
        %elseif n>0
        EC(ist,n) = C(ist,n) * E(ist,1+n) / E(ist,1+(n-1));
        %end
    end
end
%EC
%E
%Ever
%Eprime
%Eprimever
%Eprime = Eprimever;
end
