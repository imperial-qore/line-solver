%{ @file sn_get_product_form_chain_params.m
 %  @brief Extracts product-form parameters aggregated by chain
 %
 %  @author LINE Development Team
%}

%{
 % @brief Extracts product-form parameters aggregated by chain
 %
 % @details
 % This function extracts parameters from a network structure and aggregates
 % them by chain for use in product-form analysis methods. The mu parameter
 % includes extra elements beyond population |N| as required by MVALDMX.
 %
 % @par Syntax:
 % @code
 % [lambda,D,N,Z,mu,S,V] = sn_get_product_form_chain_params(sn)
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
 % <tr><td>lambda<td>Chain arrival rates
 % <tr><td>D<td>Chain service demands at queuing stations
 % <tr><td>N<td>Chain populations
 % <tr><td>Z<td>Chain think times (service demands at delay stations)
 % <tr><td>mu<td>Load-dependent service capacity scaling factors
 % <tr><td>S<td>Number of servers at queuing stations
 % <tr><td>V<td>Chain visit ratios
 % </table>
%}
function [lambda,D,N,Z,mu,S,V]= sn_get_product_form_chain_params(sn)

[lambda,~,~,~,mu,~] = sn_get_product_form_params(sn);
queueIndex = find(sn.nodetype == NodeType.Queue);
delayIndex = find(sn.nodetype == NodeType.Delay);
ignoreIndex = find(sn.nodetype == NodeType.Source | sn.nodetype == NodeType.Join);
[Dchain,~,Vchain,~,Nchain,~,~] = sn_get_demands_chain(sn);
lambda_chains = zeros(1,sn.nchains);
for c=1:sn.nchains
    lambda_chains(c) = sum(lambda(sn.inchain{c}), 'omitnan');
    %if sn.refclass(c)>0
        %D_chains(:,c) = Dchain(find(isfinite(sn.nservers)),c)/alpha(sn.refstat(c),sn.refclass(c));
        %Z_chains(:,c) = Dchain(find(isinf(sn.nservers)),c)/alpha(sn.refstat(c),sn.refclass(c));
    %else
        D_chains(:,c) = Dchain(sn.nodeToStation(queueIndex),c);
        Z_chains(:,c) = Dchain(sn.nodeToStation(delayIndex),c);
    %end
end
S = sn.nservers(sn.nodeToStation(queueIndex));
lambda = lambda_chains;
N = Nchain;
%D_chains(sn.nodeToStation(ignoreIndex),:) =[];
Vchain(sn.nodeToStation(ignoreIndex),:) =[];
D = D_chains;
Z = Z_chains;

V = Vchain;
if isempty(Z)
    Z = 0*N;
end
end
    
    
