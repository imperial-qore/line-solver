function [soujt,sts]=ctmc_simulate(Q, pi0, n)
if isempty(pi0)
   r=rand(length(Q),1); r=r/sum(r); 
end
[~,st] =  min(abs(rand-cumsum(pi0)));
F = cumsum(Q - diag(diag(Q)),2); F=F./repmat(F(:,end),1,length(F));
for i=1:n
    sts(i) = st; soujt(i) = exprnd(-1/Q(st,st));
    st =  1 + max([0,find( rand - F(st,:) > 0)]);
end

end