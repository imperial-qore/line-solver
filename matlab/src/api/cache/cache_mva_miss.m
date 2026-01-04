%{ @file cache_mva_miss.m
 %  @brief Computes miss rates using Mean Value Analysis
 %
 %  @author LINE Development Team
%}

%{
 % @brief Computes cache miss rates using Mean Value Analysis
 %
 % @details
 % This function computes global and per-item miss rates for cache models
 % using Mean Value Analysis with given item popularities and routing.
 %
 % @par Syntax:
 % @code
 % [M, Mk] = cache_mva_miss(p, m, R)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>p<td>Item popularity probabilities
 % <tr><td>m<td>Cache capacity vector
 % <tr><td>R<td>Routing probabilities
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>M<td>Global miss rate
 % <tr><td>Mk<td>Per-item miss rate
 % </table>
%}
function [M,Mk]=cache_mva_miss(p,m,R)
n=length(p);
h=length(m);
if sum(m)==0 || min(m)<0
    Mk=ones(1,n);
    M=p*Mk';
    return
end
for j=1:h
    [~,Mj]=cache_mva_miss(p,oner(m,j),R);
    for k=1:n
        w(k,j)=prod(R(1:j,k))*p(k)^j*abs(Mj(k));
    end
end
for j=1:h
    x(j) = 1/sum(abs(w(:,j)));
end
for k=1:n
    Mk(k)=1;
    for j=1:h
        Mk(k)=Mk(k)-x(j)*m(j)*w(k,j);
    end
end
Mk=abs(Mk);
M=p*Mk';
end
