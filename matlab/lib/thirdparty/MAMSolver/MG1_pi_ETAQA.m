function pi=MG1_pi_ETAQA(B,A,G,varargin)
%MG1_pi_ETATA Aggregate probability vector of a M/G/1 process 
%Using the newETAQA method [Stathopoulos, Riska, Hua, Smirni]
%
%	Return value: pi=[pi0,pi1,pi2 + pi3 + ...];
%	
%	Usage: 
%   For Continous Case:
%	pi=MG1_pi_ETAQA(B,A,G) computes the aggregate prob. vector 
%	of an M/G/1-type Markov chain with a infinitesimal 
%	generator matrix of the form
%	
%	
%		B0  B1  B2  B3  B4  ...
%		A0  A1  A2  A3  A4  ...
%	Q =	0   A0  A1  A2  A3  ...
%		0   0   A0  A1  A2  ...
%		...
%		
%	the input matrix G is the minimal nonnegative solution to the matrix 
%	equation 0 = A0 + A1 G + A2 G^2 + ... + Aamax G^amax
%	A = [A0 A1 ... Aamax] and B = [B0 B1 B2 ... Bbmax]	
%	
%	pi=MG1_pi_ETAQA([],A,G) computes the aggregate prob. vector of an M/G/1-type 
%	Markov chain with an infinitesimal generator matrix of the form
%
%		A0  A1  A2  A3  A4  ...
%		A0  A1  A2  A3  A4  ...
%	Q =	0   A0  A1  A2  A3  ...
%		0   0   A0  A1  A2  ...
%		...
%   
%   For Discrete Case:
%   pi = MG1_pi_ETAQA(B,A,G) computes the aggregatet prob. vector 
%   of an M/G/1-type Markov chain with a transition matrix P of the form
%
%	
%		B0  B1  B2  B3  B4  ...
%		A0  A1  A2  A3  A4  ...
%	P =	0   A0  A1  A2  A3  ...
%		0   0   A0  A1  A2  ...
%		...
%	
%	the input matrix G is the minimal nonnegative solution to the matrix 
%	equation G = A0 + A1 G + A2 G^2 + ... + Aamax G^amax
%	A = [A0 A1 ... Aamax] and B = [B0 B1 B2 ... Bbmax]	
%	
%	pi=MG1_pi_ETAQA([],A,G) computes the aggregate prob. vector of an M/G/1-type 
%	Markov chain with a transition matrix P of the form
%
%		A0  A1  A2  A3  A4  ...
%		A0  A1  A2  A3  A4  ...
%	P =	0   A0  A1  A2  A3  ...
%		0   0   A0  A1  A2  ...
%		...
%
%	Optional parameters:
%		Boundary: Allows solving the MG1 type Markov chain with a more
%			  general boundary:
%
%				B0  B1  B2  B3  B4  ...
%				C0  A1  A2  A3  A4  ...
%		          Q/P=	0   A0  A1  A2  A3  ...
%				0   0   A0  A1  A2  ...
%				...
%			
%			  the parameter value contains the matrix C0.
%			  The matrices C0 and B1,B2,... need not to be square.
%			  (default: C0=A0)
%		
%			    


OptionNames = [	'Boundary    '];

OptionTypes = [	'numeric'];
		

OptionValues = [];

Options=[];
for i = 1:size(OptionNames,1)
	options.(deblank(OptionNames(i,:)))=[];
end

% Default settintgs

options.Boundary=[];

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

m=size(A,1);		% number of states for the repetitive structure	   
dega=size(A,2)/m-1;	% computer amax
A0 = A(:,1:m);


if (isempty(B))
	mb=m;			% number of boundary states
	degb=dega;		% the number of B(i) - 1, bmax
	B=A;
	if (isempty(options.Boundary))
		C0 = A0;
	else
		C0 = options.Boundary;
	end
else
	mb=size(B,1);	
	if (mod(size(B,2)-mb,m) ~= 0)
		error('Matrix B has an incorrect number of columns');
	end
	if (isempty(options.Boundary))
		C0 = A0;
	else
		C0 = options.Boundary;
	end
	degb=(size(B,2)-mb)/m;	% bmax
end

if (isempty(options.Boundary) & mb ~= m)
	error('The Boundary option must be used since a dimension of B0 is not identical to A0');
end	


if (~isempty(options.Boundary))
	if (size(options.Boundary,1) ~= m || size(options.Boundary,2) ~= mb)
		error('The boundary parameter value has an incorrect dimension');
	end
end  

if (sum(B,2) - zeros(mb,1))'*ones(mb,1) > 1e-12
% if this is a transition matrix, transform it to a infinitesimal generator matrix    
% change B0 and A1 by substracting from them the identity matrix
    if ((sum(B,2) - ones(mb,1))'*ones(mb,1) < 1e-12)
        B(:,1:mb) = B(:,1:mb) - eye(mb);
        A(:,m+1:2*m) = A(:,m+1:2*m) - eye(m);    
    end
end
    



% test the drift condition
% alpha = (sum_i=-1 Ai*i)*e
% a = stochastic_root(sum_i=-1 Ai);
% drift = alpha*a

sumA = A(:,dega*m+1:end);
alpha = sum(sumA,2);
for i = dega-1:-1:1
	sumA=sumA+A(:,i*m+1:(i+1)*m);
	alpha=alpha+sum(sumA,2);
end
sumA=sumA+A(:,1:m);
a=stat(sumA);
drift=a*alpha;

if (drift >= 1)
	error('MATLAB:MG1_pi_ETAQA:NotPositiveRecurrent',...
		'The Markov chain characterized by A is not positive recurrent');
end

% begin the exact computation for boundary and repetitive aggregated probabilities
% Preparation 1: compute the S_hat
% S_hat(j) = B(j)*G^0 + B(j+1)*G + B(j+2)*G^2 + B(j+3)*G^3 + ..., j >= 1

Shat = B(:,mb + (degb-1)*m+1:end);
if degb <= 1
    % this is for the QBD case 
    % do nothing here because Shat already has a value
else
    % this is for the non QBD case
    for i = degb-1:-1:1
        temp = B(:,mb+(i-1)*m+1:mb+i*m) + Shat(:,1:m)*G;
        Shat = [temp,Shat];
    end
end
% Preparation 2: compute the S
% S(j) = A(j)*G^0 + A(j+1)*G^1 + A(j+2)*G^2 + A(j+3)*G^3 + ..., j >= 1 for the forward subscript it is 0 
S = A(:,dega*m+1:end);
if dega <=1
    error('MATLAB:MG1_pi_ETAQA:NotEnoughAblocks',...
        'The number of repeative state blocks is less than 2, this is not a irreducible Markov Chain');
else
    for i = dega-1:-1:1
        temp = A(:,i*m+1:(i+1)*m) + S(:,1:m)*G;
        S = [temp,S];
    end
end



%Preparation 3: Build Xnew, where xXnew = [1 0s]. x is the aggregated
%stationary probability vector
%Compute first block column of Xnew, Firstc = [1s]'
%Second block column of Xnew, Secondc=[B0;C0;0s]';
%Third block column of Xnew, Thirdc=[B(1)+Shat(2)*G;A(1)+S(1)*G;0s];
%Fourth block column of Xnew, Fourthc=[sum_i>=2 Bi + sum_i>=3 Shati*G;
%sum_i>=2 A(i)+sum_i>=2 S(i)*G;sum_i>=1 A(i) + sum_i>=1 S(i)*G] -1 col
Firstc=ones(mb+2*m,1);
Secondc=[B(:,1:mb);C0;zeros(m,mb)];
if (size(Shat,2) < 2*m)
	temp = zeros(mb,m);
	Shat = [Shat,temp];
end

if (size(S,2) < 2*m)
	temp = zeros(m,m);
	S = [S,temp];
end

if ((size(S,2)/m >= 2) && (size(Shat,2)/m >=2) )
    Thirdc = [B(:,mb+1:mb+m) + Shat(:,m+1:2*m)*G;A(:,m+1:2*m)+S(:,m+1:2*m)*G;zeros(m,m)];
else
    if (size(S,2)/m < 2)
        error('MATLAB:MG1_pi_ETAQA:ReducibleMarkovChain',...
            'The number of repetitive state blocks is less than 2, this is not a irreducible Markov Chain'); 
    else     
    Thirdc = [B(:,mb+1:mb+m);A(:,m+1:2*m) + S(:,m+1:2*m)*G;zeros(m,m)];
    end
end
% compute the fourth column
% Bsum = B(2) + B(3) + B(4) + ...
% Shat_sum = Shat(3) + Shat(4) + ...
% Asum = A(2) + A(3) + A(4) + ...
% Ssum = S(2) + S(2) + S(3) + ...
Bsum=zeros(mb,m);
Shat_sum=zeros(mb,m);
if (degb <= 2)
    if (degb == 1)
        Bsum = Bsum;
    else
        Bsum = Bsum + B(:,mb+1*m+1:mb+2*m);
    end
else 
    for i = 2:degb-1
        Bsum=Bsum+B(:,mb+(i-1)*m+1:mb+(i)*m);
        Shat_sum=Shat_sum+Shat(:,i*m+1:(i+1)*m);
    end
    Bsum=Bsum + B(:,mb+(degb-1)*m+1:end);
end

	


Asum = zeros(m,m);
Ssum = zeros(m,m);
if dega >= 3
    for i = 2:dega-1
        Ssum = Ssum + S(:,i*m+1:(i+1)*m);
        Asum = Asum + A(:,i*m+1:(i+1)*m);
    end
    Asum = Asum + A(:,dega*m+1:end);
else
   if dega == 2
        Asum = Asum + A(:,dega*m+1:end);
   else
       error('MATLAB:MG1_pi_ETAQA:NotEnoughA Blocks',...
           'The number of repetitive state blocks are less than 3, the Markov Chain is reducible');
   end
end


Fourthc = [Bsum+Shat_sum*G;Asum+Ssum*G;Asum+A(:,m+1:2*m)+(Ssum+S(:,m+1:2*m))*G];
Xtemp = [Secondc,Thirdc,Fourthc];

i = 0;
rankt = 0;
rankXtemp = rank(Xtemp);
while (rankt ~= rankXtemp && i <= (mb+2*m-1))
	rankt = rank([Xtemp(:,1:i),Xtemp(:,i+2:mb+2*m)]);
	i = i+1;
end



Xtemp = [Xtemp(:,1:i-1),Xtemp(:,i+1:mb+2*m)];


Xnew = [Firstc,Xtemp];





rside = [1,zeros(1,mb+2*m-1)];
pi = rside / Xnew;

end
