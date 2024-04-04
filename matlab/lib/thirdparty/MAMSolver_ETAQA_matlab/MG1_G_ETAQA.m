function G = MG1_G_ETAQA(An)
% MG1_G_ETAQA computes the G matrix for M/G/1-Type Markov Chains [Bini, Meini]
%
%
%	Discrete case:
%	G = MG1_G_ETAQA(An) computes the minimal nonnegative solution to the 
%	matrix equation G = A0 + A1 G + A2 G^2 + A3 G^3 + ... + A_max G^max,
% 	a nonnegative matrix, with (A0+A1+A2+...+A_max) irreducible and
%	stochastic
%
%	Continuous case:
%	G = MG1_G_ETAQA(An) computes the minimal nonnegative solution to the
%	matrix equation 0 = A0 + A1 G + A2 G^2 + A3 G^3 + ... + A_max G^max,
%	a nonnegative matrix, with (A0+A1+A2+...+A_max) irreducible and
%       stochastic 

r = size(An,1);
s = size(An,2);
b = s/r;
psize = r;
Gtemp = zeros(psize,psize);


% check if the matrix is stochastic

rowsum = sum(An,2);
isdiscrete = 0;
for i=1:r
	if (rowsum(i) < 1e-16)
		
	else
		isdicrete = 1;
	end
end

if (isdiscrete == 0)
	L = An(:,r+1:2*r);
	t = min(diag(L));
	if (t > 0)
		error('This is not stochastic matrix, neither continuous nor discrete! \n Please make sure every row sum up to 0 or 1');	
	end
	An = An/(-t);
	An(:,r+1:2*r) = An(:,r+1:2*r) + eye(r);

end


G = MG1_CR(An);


%while max(max(abs(G-Gtemp))) > 1e-15
%	Gtemp = G;
%	G = zeros(psize,psize);
%    E=eye(psize);
%	for i = 1:s
%		G = G + An{i}*E;
%        E=E*Gtemp;
%	end
%end




end 
