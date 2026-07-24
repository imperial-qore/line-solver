function qlen = QBD_qlen_ETAQA(B,A,B0,pi,n)
%QBD_qlen_ETAQA computes n-th moment of a Quasi-Birth-Death(QBD) process
%using the ETAQA method [Riska, Smirni]
%
%
%	Return value: qlen = E[qlen^n];
%
%	Usage:
%   	For Continuous Case:
%	qlen=QBD_qlen_ETAQA(B,A,B0,pi,n) computes the n-th moment of Q length of a 
%	QBD-type Continuous Time Markov Chain(CTMC) with an infinitesimal generator
%	matrix of the form,
%	
%           B1  B2  0   0   0   ...
%           B0  A1  A2  0   0   ...
%	 Q =  	0   A0  A1  A2  0   ...
%           0   0   A0  A1  A2  ...
%           ...
%	where B = [B1,B2] and A = [A0,A1,A2], B0 describes the backward transition from
%	second level of states to the boundary level of states. The input pi is the aggregated 
%   	probabilities computed from ETAQA method.
%   
%
%   	For Discrete Case:
%   	pi=QBD_qlen_ETAQA(B,A,B0,pi,n) computes the n-th moment of Q length  of a 
%   	QBD-type Discrete Time Markov with a transition matrix of the form,
%   
%           B1  B2  0   0   0   ...
%           B0  A1  A2  0   0   ...
%	 P =0   A0  A1  A2  0   ...
%           0   0   A0  A1  A2  ...
%           ...
%	where B = [B1,B2] and A = [A0,A1,A2], B0 describes the backward transition from
%	second level of states to the boundary level of states. The input pi is
%	the aggregated probabilities computed from ETAQA method
%	
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


if abs((sum(B,2) - ones(mb,1))'*ones(mb,1)) < 1e-12
   % if this is a transition matrix transform it to the infinitesimal
   % generator
   
   A(:,m+1:2*m) = A(:,m+1:2*m) - eye(m);
   B(:,1:mb) = B(:,1:mb) - eye(mb);
end

pi0 = pi(1:mb);
pi1 = pi(mb+1:mb+m);
pistar = pi(mb+m+1:mb+2*m);

lsleft = A(:,1:m) + A(:,m+1:2*m) + A(:,2*m+1:3*m);
lsleft = lsleft(:,1:end-1);
lsleft = [lsleft,(A(:,2*m+1:3*m) - A(:,1:m))*ones(m,1)];

r = pistar;

for k = 1:n
    bk = (-1)*(2^k*pi0*B(:,mb+1:mb+m) + 2^k*pi1*A(:,m+1:2*m) + 3^k*pi1*A(:,2*m+1:3*m));
    bkrest = zeros(1,m);
    for l = 1:k
       bkrest = bkrest + bino(k,l)*(2^l*r(k-l+1,:)*A(:,2*m+1:3*m) + r(k-l+1,:)*A(:,m+1:2*m));
    end
    bk = bk - bkrest;  
    ck = -2^k*pi1*A(:,2*m+1:3*m)*ones(m,1);
    ckrest = 0;
    for l = 1:k
        ckrest = ckrest + bino(k,l)*r(k-l+1,:)*A(:,2*m+1:3*m)*ones(m,1);
    end
    bk = bk(1:end-1);
    ck = ck - ckrest;
    rside = [bk,ck];
    rk = rside /lsleft;
    r = [r;rk];
end

qlen = sum(r(end,:),2) + sum(pi1,2);


	function b = bino(n,k)
		b = factorial(n)/(factorial(k)*factorial(n-k));
	end

end
