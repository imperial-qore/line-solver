%{
%{
 % @file pfqn_mushift.m
 % @brief Shift load-dependent service rate vector by removing first element.
%}
%}

%{
%{
 % @brief Shift load-dependent service rate vector by removing first element.
 % @fn pfqn_mushift(mu, iset)
 % @param mu Load-dependent rate matrix (MxN).
 % @param iset Set of station indices to shift.
 % @return mushifted Shifted rate matrix (Mx(N-1)).
%}
%}
function mushifted=pfqn_mushift(mu,iset)
% shifts the service rate vector
[M,N]=size(mu);

for i=iset(:)'
    for m=1:M
        if m==i
            mushifted(m,1:(N-1))=mu(m,2:N);
        else
            mushifted(m,1:(N-1))=mu(m,1:(N-1));
        end        
    end
end
end