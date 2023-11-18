function pi=MG1_pi(B,A,G,varargin)
%MG1_pi Stationary vector of a Quasi-Birth-Death Markov Chains [Neuts, Ramaswami] 
%
%   pi=MG1_pi(B,A,G) computes the stationary vector of an M/G/1-type
%   Markov chain with a transition matrix of the form
%
%               B0  B1  B2  B3  B4  ...               
%               A0  A1  A2  A3  A4  ...
%       P  =    0   A0  A1  A2  A3  ...
%               0   0   A0  A1  A2 ...
%               ...
%
%   the input matrix G is the minimal nonnegative solution to the matrix 
%   equation G = A0 + A1 G + A2 G^2 + ... + Aamax G^amax
%   A = [A0 A1 ... Aamax] and B=[B0 B1 B2 ... Bbmax]
%
%   pi=MG1_pi([],A,G) computes the stationary vector of an M/G/1-type
%   Markov chain with a transition matrix of the form
%
%               A0  A1  A2  A3  A4  ...               
%               A0  A1  A2  A3  A4  ...
%       P  =    0   A0  A1  A2  A3  ...
%               0   0   A0  A1  A2 ...
%               ...
%
%   Optional Parameters:
%   
%       MaxNumComp: Maximum number of components (default: 500)
%       Verbose: The accumulated probability mass is printed at every 
%                n steps when set to n (default:0)
%       Boundary: Allows solving the MG1 type Markov chain with a more 
%                 general boundary:
%
%                           B0  B1  B2  B3  B4  ...               
%                           C0  A1  A2  A3  A4  ...
%                   P  =    0   A0  A1  A2  A3  ...
%                           0   0   A0  A1  A2 ...
%                           ...
%                           
%                 the parameter value contains the matrix C0.
%                 The matrices C0 and B1,B2,... need not to be square.
%                 (default: C0=A0)
%       InputPiZero: The first component of the steady state vector pi
%                    can be given as input by setting this parameter value
%                    equal to its first component pi_0 (default: not given)


OptionNames=[
             'MaxNumComp ';
             'Verbose    ';
             'Boundary   ';
             'InputPiZero'];
OptionTypes=[
             'numeric';
             'numeric';
             'numeric';
             'numeric'];
 
OptionValues=[];             
             
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.MaxNumComp=500;
options.Verbose=0;
options.Boundary=[];
options.InputPiZero=[];

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

m=size(A,1);
dega=size(A,2)/m-1;

if (isempty(B))
    mb=m;
    degb=dega;
    if (~isempty(options.Boundary))
        B=A;
    end    
else
    mb=size(B,1);
    if (mod(size(B,2)-mb,m) ~= 0)
        error('Matrix B has an incorrect number of columns');
    end    
    degb=(size(B,2)-mb)/m;
end

if (isempty(options.Boundary) & mb ~= m)
    error('The Boundary option must be used as dimension of B0 is not identical to A0');
end

if (~isempty(options.Boundary))
    if (size(options.Boundary,1) ~= m || size(options.Boundary,2) ~= mb)
        error('The Boundary parameter value has an incorrect dimension');
    end    
end    
% Compute theta and beta, sum_v>=0 Av and sum_v>=k Av G^v-1 
% the last sums (for k=1,...,amax) are stored in A
sumA=A(:,dega*m+1:end);
beta=sum(sumA,2);
% beta = (A_maxd)e + (A_maxd + A_maxd-1)e + ... + (Amaxd+...+A1)e
for i=dega-1:-1:1
    sumA=sumA+A(:,i*m+1:(i+1)*m);
    A(:,i*m+1:(i+1)*m)=A(:,i*m+1:(i+1)*m)+A(:,(i+1)*m+1:(i+2)*m)*G;
    beta=beta+sum(sumA,2);
end
sumA=sumA+A(:,1:m);
theta=stat(sumA);
drift=theta*beta;

if (drift >= 1)
    error('MATLAB:MG1_pi:NotPossitiveRecurrent',...
        'The Markov chain characterized by A is not positive recurrent');
end


if (isempty(options.InputPiZero))
    % Compute g
    g=stat(G);
    
    if (isempty(B))
        % Compute pi_0
        pi0=(1-drift)*g;
    else
        % Compute sum_v>=1 Bv, sum_v>=1 (v-1) Bv e, sum_v>=k Bv G^v-1
        % the last sums (for k=1,...,bmax) are stored in B
        sumBB0=B(:,mb+(degb-1)*m+1:end);
        Bbeta=zeros(mb,1);
        for i=degb-1:-1:1
            Bbeta=Bbeta+sum(sumBB0,2);
            sumBB0=sumBB0+B(:,mb+(i-1)*m+1:mb+i*m);
            B(:,mb+(i-1)*m+1:mb+i*m)=B(:,mb+(i-1)*m+1:mb+i*m)+B(:,mb+i*m+1:mb+(i+1)*m)*G;
        end
      
        
        % Compute K, kappa
        if (isempty(options.Boundary)) %m=mb
            K=B(:,1:mb)+B(:,mb+1:mb+m)*G;
        else
            L=(eye(m)-A(:,m+1:2*m))^(-1)*options.Boundary;
            K=B(:,1:mb)+B(:,mb+1:mb+m)*L;
        end
        kappa=stat(K);

        % Compute psi1, psi2
        temp=sum((eye(m)-sumA-(ones(m,1)-beta)*g)^(-1),2);
        psi1=(eye(m)-A(:,1:m)-A(:,m+1:2*m))*temp+...
            (1-drift)^(-1)*sum(A(:,1:m),2);
        psi2=ones(mb,1)+(sumBB0-B(:,mb+1:mb+m))*temp+(1-drift)^(-1)*Bbeta;

        % Compute kappa1
        tildekappa1=psi2+B(:,mb+1:mb+m)*(eye(m)-A(:,m+1:2*m))^(-1)*psi1;

        % Compute pi_0
        pi0=(kappa*tildekappa1)^(-1)*kappa;
    end
else
    if (~isempty(B)) 
        % Compute sum_v>=1 Bv, sum_v>=1 (v-1) Bv e, sum_v>=k Bv G^v-1
        % the last sums (for k=1,...,bmax) are stored in B
        sumBB0=B(:,mb+(degb-1)*m+1:end);
        Bbeta=zeros(mb,1);
        for i=degb-1:-1:1
            Bbeta=Bbeta+sum(sumBB0,2);
            sumBB0=sumBB0+B(:,mb+(i-1)*m+1:mb+i*m);
            B(:,mb+(i-1)*m+1:mb+i*m)=B(:,mb+(i-1)*m+1:mb+i*m)+B(:,mb+i*m+1:mb+(i+1)*m)*G;
        end
    end
    pi0=options.InputPiZero;
end    

numit=1;
sumpi=sum(pi0);
pi = [];

% Start stable RAMASWAMI formula
invbarA1=(eye(m)-A(:,m+1:2*m))^(-1);
while (sumpi < 1-10^(-10) && numit < options.MaxNumComp)
    if (numit <= degb)
        if (isempty(B))%mb=m
            pi(numit,1:m)=pi0*A(:,numit*mb+1:(numit+1)*mb);
        else
            pi(numit,1:m)=pi0*B(:,mb+(numit-1)*m+1:mb+numit*m);
        end
    else
        pi(numit,1:m)=zeros(1,m);
    end
    for j=1:min(numit-1,dega-1)
        pi(numit,1:m)=pi(numit,1:m)+...
            pi(numit-j,:)*A(:,(j+1)*m+1:(j+2)*m);
    end    
    pi(numit,:)=pi(numit,:)*invbarA1;
    sumpi=sumpi+sum(pi(numit,:));
    if (~mod(numit,options.Verbose))
            fprintf('Accumulated mass of the first %d (reblocked) components: %d\n',numit,sumpi);
            drawnow;
    end
    numit=numit+1;
end    
pi=[pi0 reshape(pi',1,[])];
if (numit == options.MaxNumComp)
    warning('Maximum Number of Components %d reached',numit);
end