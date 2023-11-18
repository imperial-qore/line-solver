function pi=GIM1_pi(B,R,varargin)
%GIM1_pi Stationary vector of a GI/M/1-Type Markov Chain [Neuts] 
%
%   pi=GIM1_pi(B,R) computes the stationary vector of a GI/M/1-Type
%   Markov chain with a transition matrix of the form
%
%               B1  A0  0   0   0  ...               
%               B2  A1  A0  0   0  ...
%       P  =    B3  A2  A1  A0  0  ...
%               B4  A3  A2  A1  A0 ...
%               ...
%
%   the input matrix R is the minimal nonnegative solution to the matrix 
%   equation R = A0 + R A1 + R^2 A2 + ... + R^maxa Amaxa
%   The input matrix B equals [B1; B2; ...; Bmaxb].
%
%   Optional Parameters:
%   
%       MaxNumComp: Maximum number of components (default: 500)
%       Verbose: The accumulated probability mass is printed at every 
%                n steps when set to n (default:0)
%       Boundary: Allows solving the GI/M/1 Type MC with a more general 
%                 boundary
%
%                           B1  B0  0   0   0  ...               
%                           B2  A1  A0  0   0  ...
%                   P  =    B3  A2  A1  A0  0  ...
%                           B4  A3  A2  A1  A0 ...
%                           ...
%                 the parameter value contains the matrix [B0; A1; A2; ...; 
%                 Aamax].
%                 The matrices B0 and Bi, i > 1, need not to be square.
%                 (default: B0=A0)


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
options.MaxNumComp=500;
options.Verbose=0;

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

m=size(R,1);

temp=(eye(m)-R)^(-1);
if( max(temp<-100*eps) )
    error('MATLAB:QBD_pi:InvalidRInput',...
        'The spectral radius of R is not below 1: QBD is not pos. recurrent');
end    

if (isempty(options.Boundary))
    maxb=size(B,1)/m;
    BR=B((maxb-1)*m+1:end,:);
    for i=maxb-1:-1:1
        BR=R*BR+B((i-1)*m+1:i*m,:);
    end    
    pi=stat(BR); % compute pi_0
    pi=pi/(pi*temp*ones(m,1)); % normalize pi_0
    sumpi=sum(pi);
    numit=1;
    while (sumpi < 1-10^(-10) && numit < 1+options.MaxNumComp)
        pi(numit+1,1:m)=pi(numit,:)*R; % compute pi_(numit+1)
        numit=numit+1;
        sumpi=sumpi+sum(pi(numit,:));
        if (~mod(numit,options.Verbose))
            fprintf('Accumulated mass after %d iterations: %d\n',numit,sumpi);
            drawnow;
        end
    end   
    pi=reshape(pi',1,[]);
else
    mb=size(B,2) % number of states of boundary level
    maxbm1=(size(B,1)-mb)/m % maxb - 1
    BR1=B(mb+(maxbm1-1)*m+1:end,:);
    for i=maxbm1-1:-1:1
        BR1=R*BR1+B(mb+(i-1)*m+1:mb+i*m,:);
    end    
    maxa=(size(options.Boundary,1)-mb)/m 
    BR0=options.Boundary(mb+(maxa-1)*m+1:end,:);
    size(BR0)
    for i=maxa-1:-1:1
        BR0=R*BR0+options.Boundary(mb+(i-1)*m+1:mb+i*m,:);
    end    
    pi0=stat([B(1:mb,:) options.Boundary(1:mb,:); BR0 BR1]); % compute pi_0 and pi_1
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
    warning('Maximum Number of Components %d reached',numit-1);
end
