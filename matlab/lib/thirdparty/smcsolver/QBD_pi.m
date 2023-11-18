function pi=QBD_pi(B0,B1,R,varargin)
%QBD_pi Stationary vector of a Quasi-Birth-Death Markov Chains [Neuts] 
%
%   DISCRETE TIME CASE:
%
%   pi=QBD_pi(B0,B1,R) computes the stationary vector of a Quasi-Birth-Death
%   Markov chain with a transition matrix of the form
%
%               B1  A2  0   0   0  ...               
%               B0  A1  A2  0   0  ...
%       P  =    0   A0  A1  A2  0  ...
%               0   0   A0  A1  A2 ...
%               ...
%
%   the input matrix R is the minimal nonnegative solution to the matrix 
%   equation R = A2 + R A1 + R^2 A0
%
%   CONTINUOUS TIME CASE:
%
%   pi=QBD_pi(B0,B1,R) computes the stationary vector of a Quasi-Birth-Death
%   Markov chain with a rate matrix of the form
%
%               B1  A2  0   0   0  ...               
%               B0  A1  A2  0   0  ...
%       Q  =    0   A0  A1  A2  0  ...
%               0   0   A0  A1  A2 ...
%               ...
%
%   the input matrix R is the minimal nonnegative solution to the matrix 
%   equation 0 = A2 + R A1 + R^2 A0
%
%   Optional Parameters:
%   
%       MaxNumComp: Maximum number of components (default: 500)
%       Verbose: The accumulated probability mass is printed at every 
%                n steps when set to n (default:0)
%       Boundary: Allows solving the QBD with a more general boundary
%   
%                                 B1  B2  0   0   0  ...               
%                                 B0  A1  A2  0   0  ...
%                   P (or Q)  =   0   A0  A1  A2  0  ...
%                                 0   0   A0  A1  A2 ...
%                                 ...
%
%                 the parameter value contains the matrix [B2; A1+R*A0].
%                 The matrices B0 and B2 need not to be square.
%                 (default: B2=A2)

OptionNames=[
             'Boundary   '; 
             'MaxNumComp ';
             'Verbose    '];
OptionTypes=[
             'numeric'; 
             'numeric';
             'numeric'];
 
OptionValues=[];             
             
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Boundary=[];
options.MaxNumComp=10000;
options.Verbose=0;

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

m=size(R,1);

% Convert to discrete time problem, if needed
if (sum(diag(B1)<0)) % continuous time
    lamb=max(-diag(B1));
    if (~isempty(options.Boundary))
        mb=size(B1,1);
        lamb=max(lamb,max(-diag(options.Boundary(mb+1:end,:))));
        options.Boundary=options.Boundary/lamb;
        options.Boundary(mb+1:end,:)=options.Boundary(mb+1:end,:)+eye(m);
        B1=B1/lamb+eye(mb);
    else
        B1=B1/lamb+eye(m);
    end    
    B0=B0/lamb;
end

temp=(eye(m)-R)^(-1);
if( max(temp<-100*eps) )
    error('MATLAB:QBD_pi:InvalidRInput',...
        'The spectral radius of R is not below 1: QBD is not pos. recurrent');
end    

if (isempty(options.Boundary))
    pi=statvec(B1+R*B0); % compute pi_0
    pi=pi/(pi*temp*ones(m,1)); % normalize pi_0
    sumpi=sum(pi);
    numit=1;
    pi(2:options.MaxNumComp,:)=0;
    while (sumpi < 1-10^(-10) && numit < 1+options.MaxNumComp)
        pvec = pi(numit,:)*R;
        pi(numit+1,1:m)=pvec; % compute pi_(numit+1)
        numit=numit+1;
        sumpi=sumpi+sum(pi(numit,:));
        if (~mod(numit,options.Verbose))
            fprintf('Accumulated mass after %d iterations: %d\n',numit,sumpi);
            drawnow;
        end
    end   
%    pi=reshape(pi',1,[]);
else
    mb=size(B1,1);
    pi0=statvec([[B1; B0] options.Boundary]); % compute pi_0 and pi_1
    pi0=pi0/(pi0(1:mb)*ones(mb,1)+pi0(mb+1:end)*temp*ones(m,1)); % normalize
    pi=pi0(mb+1:end);
    pi0=pi0(1:mb);
    sumpi=sum(pi0)+sum(pi);
    numit=1;
    while (sumpi < 1-10^(-10) && numit < options.MaxNumComp)
        pi(numit+1,1:m)=pi(numit,:)*R; % compute pi_(numit+1)
        numit=numit+1;
        sumpi=sumpi+sum(pi(numit,:));
        if (~mod(numit,options.Verbose))
            fprintf('Accumulated mass after %d iterations: %d\n',numit,sumpi);
            drawnow;
        end
    end    
    pi=[pi0 reshape(pi',1,[])];
    numit=numit+1;
end    

if (numit == 1+options.MaxNumComp)
%    warning('Maximum Number of Components %d reached',numit-1);
end
