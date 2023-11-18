function [Q,U,R,T,C,X] = snDeaggregateChainResults(sn, Lchain, ST, STchain, Vchain, alpha, Qchain, Uchain, Rchain, Tchain, Cchain, Xchain)
% [Q,U,R,T,C,X] = SNDEAGGREGATECHAINRESULTS(sn, Lchain, ST, STchain, Vchain, alpha, Qchain, Uchain, Rchain, Tchain, Cchain, Xchain)
%
% Obtain class-level metrics from the chain-level ones

if isempty(ST)
    ST = 1 ./ sn.rates;
    ST(isnan(ST))=0;
end

if ~isempty(Cchain)
    error('Cchain input to snDeaggregateChainResults not yet supported');
end

Vsink = cellsum(sn.nodevisits);
Vsink = Vsink(sn.nodetype==NodeType.ID_SINK,:);

S = sn.nservers;
%mu = sn.lldscaling; 
%gamma = sn.cdscaling;

for c=1:sn.nchains
    inchain = sn.inchain{c};
    for k=inchain(:)'
        if isinf(sum(sn.njobs(inchain)))
            % open chain
            X(k) = Xchain(c) * Vsink(k); 
        else
            X(k) = Xchain(c) * alpha(sn.refstat(k),k);
        end
        for i=1:sn.nstations
            if isempty(Uchain)
                if isinf(S(i))
                    U(i,k) = ST(i,k) * (Xchain(c) * Vchain(i,c) / Vchain(sn.refstat(k),c)) * alpha(i,k);
                else
                    U(i,k) = ST(i,k) * (Xchain(c) * Vchain(i,c) / Vchain(sn.refstat(k),c)) * alpha(i,k) / S(i);
                end
            else
                if isinf(S(i))
                    U(i,k) = ST(i,k) * (Xchain(c) * Vchain(i,c) / Vchain(sn.refstat(k),c)) * alpha(i,k);
                else
                    U(i,k) = Uchain(i,c) * alpha(i,k);
                end
            end
            if Lchain(i,c) > 0
                if ~isempty(Qchain)
                    Q(i,k) = Qchain(i,c) * alpha(i,k);
                else
                    Q(i,k) = Rchain(i,c) * ST(i,k) / STchain(i,c) * Xchain(c) * Vchain(i,c) / Vchain(sn.refstat(k),c) * alpha(i,k);
                end
                T(i,k) = Tchain(i,c) * alpha(i,k);
                R(i,k) = Q(i,k) / T(i,k);
                %R(i,k) = Rchain(i,c) * ST(i,k) / STchain(i,c) * alpha(i,k) / sum(alpha(sn.refstat(k),inchain)');
            else
                T(i,k) = 0;
                R(i,k) = 0;
                Q(i,k) = 0;
            end
        end
        C(k) = sn.njobs(k) / X(k);
    end
end

Q=abs(Q);
R=abs(R);
X=abs(X);
U=abs(U);
T=abs(T);
C=abs(C);
T(~isfinite(T))=0;
U(~isfinite(U))=0;
Q(~isfinite(Q))=0;
R(~isfinite(R))=0;
X(~isfinite(X))=0;
C(~isfinite(C))=0;

end