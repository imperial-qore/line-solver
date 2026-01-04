function R = GIM1_R_ETAQA(An)
%GIM1_R determines R matrix of a GI/M/1-Type Markov Chain
%	
%	DISCRETE TIME case:
%	R = GIM1_R(An) computes the minimal nonnegative solution to the
%	matrix equation R = A0 + R A1 + R^2 A2 + R^3 A3 + ... + R^max A_max,
%	where A = [A0 A1 A2 A3 ... A_max] has m rows and m*max columns and is 
%	a nonnegative matrix, with (A0 + A1 + A2 + ... + A_max) irreducible
%	and stochastic.
%
%	

r = size(An,1);
s = size(An,2);
b = r/s;

rowsum = sum(An,2);
isdiscrete = 0;

L = An(s+1:2*s,:);
t = min(diag(L));
if ( t > 0 )
isdiscrete = 1;
end

if (isdiscrete == 0)
	L = An(s+1:2*s,:);
	t = min(diag(L));
	if (t > 0)
		error('This is not stochastic matrix, neither continuous nor discrete! \n Please make sure every row sum up to 0 or 1');

	end
	An = An/(-t);
	An(s+1:2*s,:) = An(s+1:2*s,:) + eye(s);
end

Bn = [];
for i = 1:b
	temp = An((i-1)*s+1:i*s,:);
	Bn = [Bn temp];
end


R = GIM1_R(Bn,'A','FI');

end
