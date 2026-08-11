%{ @file sn_deaggregate_chain_results.m
 %  @brief Deaggregates chain-level performance metrics to class-level metrics
 %
 %  @author LINE Development Team
%}

%{
 % @brief Deaggregates chain-level performance metrics to class-level metrics
 %
 % @details
 % This function converts chain-level performance metrics (queue lengths,
 % utilizations, response times, throughputs) to class-level metrics using
 % aggregation factors.
 %
 % @par Syntax:
 % @code
 % [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, ST, STchain, Vchain, alpha, Qchain, Uchain, Rchain, Tchain, Cchain, Xchain)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>Lchain<td>Chain demands
 % <tr><td>ST<td>Class service times
 % <tr><td>STchain<td>Chain service times
 % <tr><td>Vchain<td>Chain visits
 % <tr><td>alpha<td>Aggregation factors mapping classes to chains
 % <tr><td>Qchain<td>Chain queue lengths
 % <tr><td>Uchain<td>Chain utilizations
 % <tr><td>Rchain<td>Chain response times
 % <tr><td>Tchain<td>Chain throughputs
 % <tr><td>Cchain<td>Chain system response times
 % <tr><td>Xchain<td>Chain system throughputs
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>Class queue lengths
 % <tr><td>U<td>Class utilizations
 % <tr><td>R<td>Class response times
 % <tr><td>T<td>Class throughputs
 % <tr><td>C<td>Class system response times
 % <tr><td>X<td>Class system throughputs
 % </table>
%}
function [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, ST, STchain, Vchain, alpha, Qchain, Uchain, Rchain, Tchain, Cchain, Xchain)

if isempty(ST)
    ST = 1 ./ sn.rates;
    ST(isnan(ST))=0;
end

if ~isempty(Cchain)
    error('Cchain input to sn_deaggregate_chain_results not yet supported');
end

Vsink = cellsum(sn.nodevisits);
Vsink = Vsink(sn.nodetype==NodeType.Sink,:);

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