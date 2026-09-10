%{ @file sn_get_product_form_params.m
 %  @brief Extracts standard product-form parameters from the network structure
 %
 %  @author LINE Development Team
%}

%{
 % @brief Extracts standard product-form parameters from the network structure
 %
 % @details
 % This function extracts class-level parameters from a network structure for
 % use in product-form queueing network analysis. The mu parameter includes
 % extra elements beyond population |N| as required by MVALDMX.
 %
 % @par Syntax:
 % @code
 % [lambda,D,N,Z,mu,S,V] = sn_get_product_form_params(sn)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>lambda<td>Arrival rates for open classes
 % <tr><td>D<td>Service demands at queuing stations
 % <tr><td>N<td>Population vector
 % <tr><td>Z<td>Think times (service demands at delay stations)
 % <tr><td>mu<td>Load-dependent service capacity scaling factors
 % <tr><td>S<td>Number of servers at queuing stations
 % <tr><td>V<td>Visit ratios
 % </table>
%}
function [lambda,D,N,Z,mu,S,V]= sn_get_product_form_params(sn)

R = sn.nclasses;
N = sn.njobs;
queueIndices = find(sn.nodetype == NodeType.Queue);
delayIndices = find(sn.nodetype == NodeType.Delay);
sourceIndex = find(sn.nodetype == NodeType.Source);
Mq = length(queueIndices); % number of queues
Mz = length(delayIndices); % number of delays
lambda = zeros(1,R);
S = sn.nservers(sn.nodeToStation(queueIndices));
for r=1:R
    if isinf(N(r))
        lambda(r) = sn.rates(sn.nodeToStation(sourceIndex),r);
    end
end

D = zeros(Mq,R);
Nct = sum(N(isfinite(N)));
mu = ones(Mq, ceil(Nct)+max(S(isfinite(S))));
for ist=1:Mq
    for r=1:R
        c = find(sn.chains(:,r),1);
        if sn.refclass(c)>0
            D(ist,r) = sn.visits{c}(sn.nodeToStateful(queueIndices(ist)),r) / sn.rates(sn.nodeToStation(queueIndices(ist)),r) / sn.visits{c}(sn.stationToStateful(sn.refstat(r)),sn.refclass(c));
        else
            D(ist,r) = sn.visits{c}(sn.nodeToStateful(queueIndices(ist)),r) / sn.rates(sn.nodeToStation(queueIndices(ist)),r);
        end
    end
    mu(ist,1:size(mu,2)) = min(1:size(mu,2), sn.nservers(sn.nodeToStation(queueIndices(ist))));
end
Z = zeros(max(1,Mz),R);
for ist=1:Mz
    for r=1:R
        c = find(sn.chains(:,r),1);
        if sn.refclass(c)>0
            Z(ist,r) = sn.visits{c}(sn.nodeToStateful(delayIndices(ist)),r) / sn.rates(sn.nodeToStation(delayIndices(ist)),r) / sn.visits{c}(sn.stationToStateful(sn.refstat(r)),sn.refclass(c));
        else
            Z(ist,r) = sn.visits{c}(sn.nodeToStateful(delayIndices(ist)),r) / sn.rates(sn.nodeToStation(delayIndices(ist)),r);
        end
    end
end
D(isnan(D))=0;
Z(isnan(Z))=0;
V=cellsum(sn.visits);
end
