% res = LevelDependentFluidStationaryDistr (masses, iniF, KF, cloF, iniB, KB, cloB, what, points)
% calculates the startionary distribution of first order and second order
% level dependent fluid models
%
% * masses: list of masses (K+1 vectors, if there are K levels)
% * iniF, KF, cloF: initial vector, matrix exponent and closing vector for
%   the forward direction for each level, a list of size K
% * iniB, KB, cloB: initial vector, matrix exponent and closing vector for
%   the backward direction for each level, a list of size K
% * what: string, can be 'pdf', 'cdf' and 'cdfm'
%     - 'pdf'  returns the (state-dependent) density for each point
%     - 'pdfd' returns the derivative of the density
%     - 'cdf'  returns the distribution function, P(X<p)
%     - 'cdfm' returns the distribution function, P(X<=p)
%
function res = LevelDependentFluidStationaryDistr (masses, iniF, KF, cloF, iniB, KB, cloB, T, what, points)

    K = length(T);
    T = [0 T];
    N = length(masses{1});

    function KAi = integExp(KA,L)
        l=CRPSolve(KA);
        r=CRPSolve(KA')';  
        l=l/(l*r);
        KAi = inv(-(KA-r*l)) *(eye(size(KA,1))-expm((KA-r*l)*L)) + r*l*(L+exp(-L)-1);       
    end

    function [KAi,KBi] = integExp2(KA,KB,L)
        if min(abs(eig(KA))) > min(abs(eig(KB)))
            KAi = inv(-KA)*(eye(size(KA,1))-expm(KA*L));
            KBi = integExp(KB, L);
        else
            KAi = integExp(KA, L);
            KBi = inv(-KB)*(eye(size(KB,1))-expm(KB*L));
        end
    end    
    
    cummulate = strcmp(what,'cdf') || strcmp(what,'cdfm');
    res = [];
    for p=points
        pres = zeros(1,N);
        k=0;
        while k<K && p>=T(k+1)
            if cummulate
                if k>0
                    [sumKF, sumKB] = integExp2(KF{k}, KB{k}, T(k+1)-T(k));
                    val = iniF{k}*sumKF*cloF{k} + iniB{k}*sumKB*cloB{k};
                    pres = pres + val;
                end
                if p>T(k+1) || strcmp(what,'cdfm')
                    pres = pres + masses{k+1};
                end
            end
            k = k + 1;
        end
        if k==K && p==T(k+1) && strcmp(what,'cdfm')
            pres = pres + masses{k+1};       
        end
        prem = p - T(k);
        Tk = T(k+1)-T(k);
        if strcmp(what,'pdf')
            pres = iniF{k}*expm(KF{k}*prem)*cloF{k} + iniB{k}*expm(KB{k}*(Tk-prem))*cloB{k};
        elseif strcmp(what,'pdfd')
            pres = iniF{k}*KF{k}*expm(KF{k}*prem)*cloF{k} - iniB{k}*KB{k}*expm(KB{k}*(Tk-prem))*cloB{k};
        elseif strcmp(what,'cdf') || strcmp(what,'cdfm')
            [sumKF, sumKB] = integExp2(KF{k}, KB{k}, prem);
            pres = pres + iniF{k}*sumKF*cloF{k} + iniB{k}*expm(KB{k}*(Tk-prem))*sumKB*cloB{k};
        end
        res = [res; pres];
    end
end

