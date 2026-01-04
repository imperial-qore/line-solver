function qlen=MG1_qlen_ETAQA(B,A,pi,n,varargin)
%MG1_qlen_ETAQA returns the n-th moment of queue length of a
%M/G/1 process using the ETAQA method [Riska, Smirni]
%	
%
%	Usage:
%	qlen=MG1_qlen_ETAQA(B,A,pi,n) computes the n-th
%	moment of an M/G/1-type Markov chain with an 
%	infinitesimal generator matrix of the form
%
%			B0  B1  B2  B3  B4  ...
%			A0  A1  A2  A3  A4  ...
%	   Q =	0   A0  A1  A2  A3  ...
%			0   0   A0  A1  A3  ...
%			...
%
%	the input matrix A = [A0 A1 ... Aamax] and B = [B0 B1 
%	B2 ... Bbmax], pi is the aggregated probability computed from ETAQA
%	method
%
%	pi=MG1_qlen_ETAQA([],A,pi,n) computes the stationary 
%	vector of an M/G/1-type Markov chain with an 
%	infinitesimal generator matrix of the form
%
%			A0  A1  A2  A3  A4  ...
%			A0  A1  A2  A3  A4  ...
%	  Q =	0   A0  A1  A2  A3  ...
%			0   0   A0  A1  A2  ...
%			...
%
%	Optional Parameters:
%
%		Boundary: Allows solving the MG1 type Markov chain
%			    with a more general boundary:
%					B0  B1  B2  B3  B4  ...
%					C0  A1  A2  A3  A4  ...
%				Q=	0   A0  A1  A2  A3  ...
%					0   0   A0  A1  A2  ...
%					...
%			    the parameter value contains the matrix C0.
%			    The matrices C0 and B1,B2,... need not to be 
%			    square. (default: C0=A0)
%		

OptionNames = ['Boundary  '];
OptionTypes = ['numeric'];
OptionValues = [];

options=[];
for i=1:size(OptionNames,1)
	options.(deblank(OptionNames(i,:)))=[];
end
	
% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% Sanity Checks
m=size(A,1);	% number of states for the repetitive structure
dega=size(A,2)/m-1; % amax, # of A blocks - 1 

if (isempty(B))
	mb=m;		% number boudary states
	degb=dega;	% bmax
	B=A;
else
	mb=size(B,1);
	if (mod(size(B,2)-mb,m) ~= 0)
		error('Matrix B has an incorrect number of columns');
	end
	degb = (size(B,2)-mb)/m; % bmax
end

if (isempty(options.Boundary))
	C0 = A(:,1:m);
else
	C0 = options.Boundary;
end

if (isempty(options.Boundary) & mb ~= m)
	error('The Boundary option must be used since the column size of B0 is not identical to A0');
end

if (~isempty(options.Boundary))
	if (size(options.Boundary,1) ~= m | size(options.Boundary,2) ~= mb)
		error('The boundary parameter value has an incorrect dimension');
	end
end

if ( abs(sum(pi,2)-1) > 1e-10 ) 
	error('The input probability vector does not sum up to 1');
end

if ( mod(size(pi,2)-mb,m) ~=0 )
	error('Matlab:MG1_qlen_ETAQA:IncorrectDimension',...
	'The probability vector has an incorrect number of columns');
end
pi0 = pi(:,1:mb);
pi1 = pi(:,mb+1:mb+m);
pistar = pi(:,mb+m+1:mb+2*m);


% Provide the linear system r[k]*[(A0+A1+A2+...),(A(2)+2A(3)+3A(4)+...)*e]=[b[k],c[k]];
% b[k] = - fhat(k) - f(k) - sum_l=(1 to k) binomial(k,l)*r[k-l]*(A1+sum_j=1 (j+1)^l*A(j+1))  
% c[k] = - fchat(k) - fc(k) - sum_l=1,k binomial(k,l)*r[k-l]*sum_j=1 j^l F(0,j)*e
Asum = A(:,1:m);
for i = 1:dega
	Asum = Asum + A(:,i*m+1:(i+1)*m);
end
F11 = A(:,2*m+1:3*m);
for i = 3:dega
	F11=F11+(i-1)*A(:,i*m+1:(i+1)*m);
end
Asum = Asum(:,1:end-1);
lsleft = [Asum,(F11-A(:,1:m))*ones(m,1)];



% Fhat0(j) = sum_l=j  B(l), j>=1
Fhat0j=B(:,mb+(degb-1)*m+1:end);
for j = degb-1:-1:1
	temp = B(:,mb+(j-1)*m+1:mb+j*m) + Fhat0j(:,1:m);
	Fhat0j = [temp,Fhat0j];
end

% F0(j) = sum_l=j A(l+1), j>=1
F0j = A(:,dega*m+1:end);
for j = dega-1:-1:2
	temp = A(:,j*m+1:(j+1)*m) + F0j(:,1:m);
	F0j = [temp,F0j];
end

r = pistar;

% preparation for compute the repetitive part in the righthand side of the equation
% frestsaver and fckrestsaver are the components reusable
% frestsaver(l) = sum_j=1 (j+1)^l*A(j+1), l = 1:n
% fckrestsaver(l) = sum_j=1 j^l*F0j(j), l = 1:n
frestsaver = {};
fcrestsaver = {};
for l = 1:n
	temp = zeros(m,m);
	for j = 2:dega
		temp = temp + j^l*A(:,j*m+1:(j+1)*m);
	end
	frestsaver{end+1} = temp;
	temp = zeros(m,m);
	for j = 1:dega
		if (j <= size(F0j,2)/m) 
			temp = temp + j^l*F0j(:,(j-1)*m+1:j*m);	
		end
	end
	fcrestsaver{end+1} = temp;
end

for k = 1:n
	% fhat(k) = pi0 * (sum_j=1 (j+1)^k*B(j))
	fhatk = zeros(mb,m);
	for j = 1:degb
		fhatk = fhatk + (j+1)^k*B(:,mb+(j-1)*m+1:mb+j*m);
	end
	fhatk = pi0*fhatk;
	% f(k) = pi1* (2^k*A(1) + sum_j=1 (j+2)^k*A(j+1))
	fk = zeros(m,m);
	for j = 2:dega
		fk = fk + (j+1)^k*A(:,j*m+1:(j+1)*m);
    end
    fk = 2^k*A(:,m+1:2*m) + fk;
	fk = pi1*fk;
	% frest = sum_l=(1 to k) bino(k,l)r[k-l]*(A(1)+sum_j=1 (j+1)^l A(j+1)) 
	frest = zeros(1,m);
	for l = 1:1:k
		frest = frest + bino(k,l)*r(k-l+1,:)*(A(:,m+1:2*m) + frestsaver{l});
	end
	% bktemp = -fhatk - fk - frest
	% bk = bktemp - last column	
	bk = (-1)*fhatk + (-1)*fk + (-1)*frest;
	bk = bk(:,1:end-1);
	% fchatk = pi0*(sum_j=2 j^k Fhat0j(j)*e)
	% fck = pi1*(sum_j=1 (j+1)^k*F0j(j)*e)
	% fcrest = sum_l=1 (bino(k,l)r[k-l]*sum_j=1 j^l*F0j(j)*e)
	fchatk= zeros(mb,m);
	for j=2:degb
		fchatk = fchatk + (j^k)*Fhat0j(:,(j-1)*m+1:j*m);
	end
	fchatk = pi0*fchatk*ones(m,1);
	fck = zeros(m,m);	
	for j=1:(dega-1)
		fck = fck + (j+1)^k*F0j(:,(j-1)*m+1:j*m); 
	end
	fck = pi1*fck*ones(m,1);
	fcrest = zeros(1,m);
	for l=1:k
		fcrest = fcrest + bino(k,l)*r(k-l+1,:)*fcrestsaver{l};
	end
	fcrest = fcrest*ones(m,1);
	ck = - fchatk - fck - fcrest;
	rside = [bk, ck];
	rk = rside / lsleft;
	r = [r;rk];
end



	function b = bino(n,k)
		b = factorial(n)/(factorial(k)*factorial(n-k));
	end

%Compute Q length
qlen = sum(r(end,:),2) + sum(pi1,2);



end
