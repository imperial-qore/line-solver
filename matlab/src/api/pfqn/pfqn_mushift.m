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