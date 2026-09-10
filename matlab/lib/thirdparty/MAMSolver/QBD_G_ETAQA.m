function G=QBD_G_ETAQA(An)
% QBD_G_ETAQA computes the G matrix for computation of truncated probabilities by ETAQA
% 
%	DISCRETE TIME case:
%	G = QBD_G_ETAQA computes the minimal nonnegative solution to the
%   	matrix equation G = A0 + A1 G + A2 G^2, where Ai are square nonnegative matrices, 
%       with An = [A0,A1,A2] and A0+A1+A2 as irreducible and stochastic
%
%       CONTINUOUS TIME case:
%       G = QBD_G_ETAQA computes the minimal nonnegative solution to the
%   	matrix equation 0 = A0 + A1 G + A2 G^2, where Ai are square nonnegative matrices, 
%       with An = [A0,A1,A2] and A0+A1+A2 as irreducible and stochastic
%
%      

t = size(An,1);
if (size(An,2) ~= 3*t)
	error('The size of An is not correct. An = [A0 A1 A2] with Ai as square matrices');
end
A0 = An(:,1:t);
A1 = An(:,t+1:2*t);
A2 = An(:,2*t+1:3*t);

G = QBD_CR(A0,A1,A2);

end

