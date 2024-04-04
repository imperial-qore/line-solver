function pi=QBD_pi_ETAQA(B,A,B0,G)
%QBD_pi_ETAQA Aggregate Stationary vector of a Quasi-Birth-Death(QBD) process
%using the newETAQA method [Stathopoulos, Riska, Hua, Smirni]
%
%
%	Return value: pi=[pi0,pi1,pi2+pi3+...];
%
%	Usage:
%   For Continuous Case:
%	pi=QBD_pi_ETAQA(B,A,B0,G) computes the aggregated probability vector of a 
%	QBD-type Continuous Time Markov Chain(CTMC) with an infinitesimal generator
%	matrix of the form,
%	
%           B1  B2  0   0   0   ...
%           B0  A1  A2  0   0   ...
%	 Q =  	0   A0  A1  A2  0   ...
%           0   0   A0  A1  A2  ...
%           ...
%	where B = [B1,B2] and A = [A0,A1,A2], B0 describes the backward transition from
%	second level of states to the boundary level of states. The input matrix R is 
%   the minimal nonnegative solution to the matrix equation 0 = A2 + R A1 + R^2 A0 
%   for the continuous case
%
%   For Discrete Case:
%   pi=QBD_pi_ETAQA(B,A,B0,G) computes the aggregated probability vector of a 
%   QBD-type Discrete Time Markov with a transition matrix of the form,
%   
%           B1  B2  0   0   0   ...
%           B0  A1  A2  0   0   ...
%	 P =  	0   A0  A1  A2  0   ...
%           0   0   A0  A1  A2  ...
%           ...
%	where B = [B1,B2] and A = [A0,A1,A2], B0 describes the backward transition from
%	second level of states to the boundary level of states. The input matrix R is 
%   the minimal nonnegative solution to the matrix equation R = A2 + R A1 + R^2 A0 
%   for the continuous case
%   

mb = size(B,1);
m = size(A,1);
testm = size(B0,1);
testmb = size(B0,2);
if (testm ~= m)
	error('Matrix B0 has an incorrect number of rows');
end
if (testmb ~= mb)
	error('Matrix B0 has an incorrect number of columns');
end
if (mod(size(B,2)-mb,m) ~= 0 || (size(B,2)-mb)/m ~=1 )
	error('Matrix B has an incorrect number of columns');
end 
if (mod(size(A, 2),m) ~= 0 || size(A,2)/m ~= 3)
	error('Matrix A has an incorret number of columns');
end
B1 = B(:,1:mb);
B2 = B(:,mb+1:mb+m);
A0 = A(:,1:m);
A1 = A(:,m+1:2*m);
A2 = A(:,2*m+1:3*m);


if abs((sum(B,2) - ones(mb,1))'*ones(mb,1)) < 1e-12
   % if this is a transition matrix transform it to the infinitesimal
   % generator
   
   A1 = A1 - eye(m);
   B1 = B1 - eye(mb);
end

% Construct X so that pi*X = [1,0s]
Firstc = [B1;B0;zeros(m,mb)];
Secondc = [B2;A1+A2*G;zeros(m,m)];
Thirdc = [zeros(mb,m);A2;A1+A2+A2*G];
Xtemp = [Firstc,Secondc,Thirdc];
rankXtemp = rank(Xtemp);
i = 0;
rankt = 0;
while ( rankt ~= rankXtemp && i <= (mb+2*m-1))
	rankt = rank([Xtemp(:,1:i),Xtemp(:,i+2:mb+2*m)]); 
	i = i+1;	
end

X = [ones(mb+2*m,1),Xtemp(:,1:i-1),Xtemp(:,i+1:mb+2*m)];

rside = [1,zeros(1,mb+2*m-1)];

pi = rside / X;


end


