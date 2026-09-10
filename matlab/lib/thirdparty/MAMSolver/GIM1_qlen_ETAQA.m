function qlen=GIM1_qlen_ETAQA(B,A,R,pi, n,varargin)
%GIM1_qlen_ETAQA returns the n-th moment of queue length of a
%GIM1/M/1 process using the ETAQA method [Riska, Smirni]
%	
%
%	Usage:
%	qlen=GIM1_qlen_ETAQA(B,A,pi,R,n) computes the n-th
%	moment of an GI/M/1-type Markov chain with an 
%	infinitesimal generator matrix of the form
%
%	
%		B1  A0  0   0   0   ...
%		B2  A1  A0  0   0   ...
%	 Q =B3  A2  A1  A0  0   ...
%		B4  A3  A2  A1  A0  ...
%		...
%
%	the input matrix R is the minimal nonnegative solution to the
%	matrix equation 0 = A0 + R A1 + R^2 A2 + ... + R^maxa Amaxa
%	The input matrix B equals [B1; B2; ...; Bmaxb], the input 
%	matrix A equals [A0;A1;A2; ...; Amaxa]
%
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

OptionNames = ['Boundary  '];
OptionTypes = ['numeric'];
OptionValues = [];

options=[];
for i=1:size(OptionNames,1)
	options.(deblank(OptionNames(i,:)))=[];
end
	
% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

m=size(R,1);
mb=size(B,2);

if (mod(size(B,1)-mb,m) ~= 0)
	error('MATLAB:GIM1_qlen_ETAQA: input matrix B has incorrect number of rows');
else
	degb = (size(B,1)-mb)/m;
end

if (mod(size(A,1),m) ~= 0)
	error('MATLAB:GIM1_qlen_ETAQA: input matrix A has incorrect number of rows');
else
	dega = (size(A,1))/m-1;
end

temp=(eye(m)-R)^(-1);
if( max(temp<-100*eps))
	error('MATLAB:GIM1_qlen_ETAQA:InvalidRInput',...
		'The spectral radius of R is not below 1: GIM1 is not pos. recurrent');
end

if (~isempty(options.Boundary))
	B0 = options.Boundary;
	testmb = size(B0,1);
	testm = size(B0,2);
	if (testmb ~= mb)
		error('MATLAB:GIM1_qlen_ETAQA:input options.Boundary has incorrect number of rows');
	end 
	if (testm ~= m)
		error('MATLAB:GIM1_qlen_ETAQA:input options.Boundary has incorrect number of columns');
	end
else
	B0 = A(1:m,:);
end

pi0 = pi(1:mb);
pi1 = pi(mb+1:mb+m);
pistar = pi(mb+m+1:mb+2*m);

if (n==0)
	qlen = 1;
	return;
end

test = B(1:mb,:)+B0;

if abs((sum(test,2) - ones(mb,1))'*ones(mb,1)) < 1e-10
    % transform from DTMC to CTMC
        B(1:mb,:) = B(1:mb,:) - eye(mb);
        A(m+1:2*m,:) = A(m+1:2*m,:) - eye(m);

end

k = 1;
Rpower = eye(m);
lsum = A(1:m,:) + A(m+1:2*m,:);
for i = 2:dega
    lsum = lsum + Rpower*A(i*m+1:(i+1)*m,:);
    Rpower = R*Rpower;
end




leftr = zeros(m,mb);
leftr_part1 = zeros(m,m);
leftr_part2 = zeros(m,m);
if degb >= 2 && dega >= 2
    Rpower = R;
    leftr_part1 = A(3);
    for i = 1:dega-2
        leftr_part1 = leftr_part1 + Rpower*A((i+2)*m+1:(i+3)*m,:);
        leftr_part2 = leftr_part2 + i*Rpower*A((i+2)*m+1:(i+3)*m,:);
        Rpower = R*Rpower;
    end
    
    leftr = leftr_part1*ones(m,1) + leftr_part2*ones(m,1) - A(1:m,:)*ones(m,1);

else
    if (degb == 1 && dega ~= 1)
        Rpower = eye(m); 
        for i = 2:dega
            leftr = leftr + (i-1)*Rpower*A(i*m+1:(i+1)*m,:);
            Rpower = R*Rpower;
        end
        leftr = leftr*ones(m,1) - A(1:m,:)*ones(m,1);
    else
        error('MATLAB:GIM1_qlen_ETAQA: the number of A blocks is not enough, this is a reducible Markov Chain!');
    end
end




% lsum is one matrix and leftr is one column 
lsleft = [lsum(:,1:end-1),leftr];



r = pistar;

% newly added here, calculate the unchanged part in the+ factor on the right
% it is -sum_j=2 to inf {(j-1)*R^(j-1)*Bhat(j)}
tempsum = zeros(m,mb);
Rpower = R;
for i=1:degb-1
	tempsum = tempsum + i*Rpower*getVerticalBlock(B,m,(i+2));
	Rpower = R*Rpower;
end
 rsidepart22 = pi1*(-tempsum*ones(mb,1));  
% compute the reusable part of the right hand side of this equation 




% have to check n to see if it is greater than 1
% get the reusable array, every column corresponds to a vector, starting from j=1
% F - sum_j=1 to inf {(sum_m=1 to j{m^l}) R^(j)*B^(j+1)}
rsidepart2reuse = [];
tempsum = zeros(m,m);
Rpower = R;
for i=1:degb-2
	tempsum = tempsum + i*Rpower*getVerticalBlock(A,m,i+3);
end





for k = 1:n
   % bk is the first part of the right side of the equation
   bk = (-1)*(pi0*2^k*B0 + pi1*(2^k*A(m+1:2*m,:) + 3^k*A(1:m,:)));
   
   for l = 1:k
       bk = bk + (-1)*bino(k,l)*r(k-l+1,:)*(A(m+1:2*m,:) + 2^l*A(1:m,:));
   end
 
   bk = bk(:,1:end-1);

   % now compute the second part of the equations system
   tempsum = zeros(m,m);
	Rpower = R;
    
	for i=1:dega-2
		temp = 0;
		for z=1:i
		   temp = temp+z^(k);	
        end

		tempsum = tempsum + temp*Rpower*getVerticalBlock(A,m,i+3);
        Rpower = Rpower*R;
		
	end
	rsidepart2reuse(:,end+1)=(A(1:m,:)-tempsum)*ones(m,1);
    
   ck = pi1*2^k*A(1:m,:)*ones(m,1);
   tempsum2 = zeros(m,mb);
   Rpower = R;
   for i=2:degb
       temp = 0;
       for z=2:i
           temp =temp + z^k;
       end
       tempsum2 = tempsum2 + temp*Rpower*B(mb + (i-1)*m+1:mb + (i)*m,:); 
       Rpower = R*Rpower;
   end
   ckresttemp = pi1*tempsum2*ones(mb,1);
    ck = ck - ckresttemp;
   ckrest = 0;

   for l = 1:k
       
       ckrest = ckrest + bino(k,l) * r(k-l+1,:)*(rsidepart2reuse(:,l));
   end
   ck = ck + ckrest;
  
   rside = [bk,ck];
   rk = rside / lsleft;
   r = [r;rk];
end

qlen = sum(r(end,:),2) + sum(pi1,2);



    function b = bino(n,k)
		b = factorial(n)/(factorial(k)*factorial(n-k));
    end
    
   function m = getVerticalBlock(Mn,blockrow,i)
                
 
				m = Mn(blockrow*(i-1)+1:blockrow*i,:);
    end

end
