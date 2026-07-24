function pi=GIM1_pi_ETAQA(B,A,R,varargin)
%GIM1_pi Aggregate Probability Vector of a GI/M/1 process
%implemented using ETAQA method [Riska, Smirni]
%
%	Return value: pi=[pi0,pi1,pi2 + pi3 + ...];
%   For Continous Case:
%	pi=GIM1_pi(B,A,R) computes the aggregate probability vetor
%	of a GI/M/1-Type Continuous Time Markov Chain(CTMC) with
%	an infinitesimal generator matrix of the form
%	
%		B1  A0  0   0   0   ...
%		B2  A1  A0  0   0   ...
%	 Q =	B3  A2  A1  A0  0   ...
%		B4  A3  A2  A1  A0  ...
%		...
%
%	the input matrix R is the minimal nonnegative solution to the
%	matrix equation 0 = A0 + R A1 + R^2 A2 + ... + R^maxa Amaxa
%	The input matrix B equals [B1; B2; ...; Bmaxb], the input 
%	matrix A equals [A0;A1;A2; ...; Amaxa]
%
%   For discrete Case:
%	pi=GIM1_pi(B,A,R) computes the aggregate probability vetor
%	of a GI/M/1-Type Discrete Time Markov Chain(DTMC) with
%	an infinitesimal generator matrix of the form
%	
%		B1  A0  0   0   0   ...
%		B2  A1  A0  0   0   ...
%	 P =	B3  A2  A1  A0  0   ...
%		B4  A3  A2  A1  A0  ...
%		...
%
%	the input matrix R is the minimal nonnegative solution to the
%	matrix equation R = A0 + R A1 + R^2 A2 + ... + R^maxa Amaxa
%	The input matrix B equals [B1; B2; ...; Bmaxb], the input 
%	matrix A equals [A0;A1;A2; ...; Amaxa]
%
%	Optional Parameters:
%
%		Boundary: B0, which allows solving the GI/M/1 Type CTMC
%			    with a more general boundary states block
%					
%					B1  B0  0   0   0   ...
%					B2  A1  A0  0   0   ...
%				Q =	B3  A2  A1  A0  0   ...
%					B4  A3  A2  A1  A0  ...
%					...
%			     The matrices B0 need not to be
%		 	     square. (default: B0=A0)
%

OptionNames = [
		    'Boundary   '];

OptionTypes = [
		     'numeric'];
OptionValues = [];
options = [];
for i = 1:size(OptionNames, 1)
	options.(deblank(OptionNames(i,:)))=[];
end

options.Boundary=[];

options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

m=size(R,1);
mb=size(B,2);

if (mod(size(B,1)-mb,m) ~= 0)
	error('MATLAB:GIM1_pi_ETAQA: input matrix B has incorrect number of rows');
else
	degb = (size(B,1)-mb)/m;
end

if (mod(size(A,1),m) ~= 0)
	error('MATLAB:GIM1_pi_ETAQA: input matrix A has incorrect number of rows');
else
	dega = (size(A,1))/m-1;
end

temp=(eye(m)-R)^(-1);
if( max(temp<-100*eps))
	error('MATLAB:GIM1_pi_ETAQA:InvalidRInput',...
		'The spectral radius of R is not below 1: GIM1 is not pos. recurrent');
end





if (~isempty(options.Boundary))
	B0 = options.Boundary;
	testmb = size(B0,1);
	testm = size(B0,2);
	if (testmb ~= mb)
		error('MATLAB:GIM1_pi_ETAQA:input options.Boundary has incorrect number of rows');
	end 
	if (testm ~= m)
		error('MATLAB:GIM1_pi_ETAQA:input options.Boundary has incorrect number of columns');
	end
else
	B0 = A(1:m,:);
end

test = B(1:mb,:)+B0;

if abs((sum(test,2) - ones(mb,1))'*ones(mb,1)) < 1e-10
    % transform from DTMC to CTMC
        B(1:mb,:) = B(1:mb,:) - eye(mb);
        A(m+1:2*m,:) = A(m+1:2*m,:) - eye(m);

end





% compute X so that pi*X = [1,0s]

Firstc = ones(mb+2*m,1);
% compute sum_i=2 R^(i-2)*(I-R)*B(i)
temp = eye(m) - R;
tempsum = zeros(m,mb);
for i = 2:degb
	tempsum = tempsum + temp*B(mb+(i-1)*m+1:mb+i*m,:);
	temp = R*temp;
end
% compute sum_i=1 R^(i-1)*(I-R)*A(i+1);
Secondc = [B(1:mb,:);B(mb+1:mb+m,:);tempsum];
temp = eye(m)-R;
tempsum = zeros(m,m);
for i = 2:dega
	tempsum = tempsum + temp*A(m+(i-1)*m+1:m+i*m,:);
	temp = R*temp;
end
Thirdc = [B0;A(m+1:m+m,:);tempsum];
%compute sum_i=1 R^i A(i+1)
temp = R;
tempsum = zeros(m,m);
for i=2:dega
	tempsum = tempsum + temp*A(m+(i-1)*m+1:m+i*m,:);
	temp = R*temp; 
end 
Fourthc = [zeros(mb,m);A(1:m,:);(A(1:m,:)+A(m+1:m+m,:)+tempsum)];
Fourthc= Fourthc(:,1:end-1);

X = [Firstc, Secondc, Thirdc, Fourthc];
rside = [1,zeros(1,mb+2*m-1)];
pi = rside /X;
	

end

