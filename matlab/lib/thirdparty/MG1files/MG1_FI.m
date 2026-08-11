function G=MG1_FI(A,varargin)
%MG1_FI Functional Iterations for M/G/1-Type Markov Chains [Neuts] 
%
%   G=MG1_FI(A) computes the minimal nonnegative solution to the 
%   matrix equation G = A0 + A1 G + A2 G^2 + A3 G^3 + ... + A_max G^max, 
%   where A = [A0 A1 A2 A3 ... A_max] has m rows and m*max columns and is
%   a nonnegative matrix, with (A0+A1+A2+...+A_max) irreducible and 
%   stochastic
%
%   G=MG1_FI(A,'NonZeroBlocks',vec) computes the minimal nonnegative 
%   solution to the matrix equation G = A0 + Av1 G^v1 + Av2 G^v2 + ...
%   Av3 G^v3 + ... + A_vmax G^vmax, where A = [A0 Av1 Av2 Av3 ... A_vmax] 
%   has m rows and m*max columns and is a nonnegative matrix, with 
%   (A0+Av1+Av2+...+A_vmax) irreducible and stochastic and vec = [v1 v2
%   v3 v4 ... vmax] is an increasing list of integers with v1 > 0
%
%   Optional Parameters:
%   
%       MaxNumIt: Maximum number of iterations (default: 10000)
%       Mode: 'Traditional': G(n+1) = (I-A1)^(-1) * ... 
%                 (A0 + A2*G(n)^2 + A3*G(n)^3 + ... + A_max*G(n)^max)
%             'Natural': G(n+1) = A0 + A1*G(n) + A2*G(n)^2+ ...
%                 A3*G(n)^3 + ... + A_max*G(n)^max)
%             'U-Based': G(n+1) = ...
%                 (I-A1-A2*G(n)-A3*G(n)^2-...-A_max*G(n)^(max-1))^(-1)*A0
%             'Shift<Mode>': where <Mode> is Traditional, Natural or
%                 U-Based uses the Shift Technique
%             (default:'U-based') 
%       Verbose: When set to k, the residual error is printed every 
%                k steps (default:0)
%       ShiftType: 'one'
%                  'tau'
%                  'dbl'
%       StartValue: Starting value for iteration (default: 0)
%       NonZeroBlocks: Reduces the computation time when several of the 
%                      Ai, i>0, matrices are equal to zero. A vector vec 
%                      passed as an argument specifies the indices of the 
%                      nonzero blocks (default: vec=[1 2 ... max])
%                      When set the possible Shift in the Mode is ignored  


OptionNames=['Mode         '; 
             'MaxNumIt     ';
             'Verbose      ';
             'ShiftType    ';
             'StartValue   ';
             'NonZeroBlocks'];
OptionTypes=['char   '; 
             'numeric';
             'numeric';
             'char   ';
             'numeric';
             'numeric'];
OptionValues{1}=['Traditional      '; 
                 'Natural          ';
                 'U-Based          ';
                 'ShiftTraditional ';
                 'ShiftNatural     ';
                 'ShiftU-Based     '];
OptionValues{4}=['one';
                 'tau';
                 'dbl'];
                 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Mode='U-Based';
options.MaxNumIt=10000;
options.ShiftType='one';
options.Verbose=0;
m=size(A,1);
options.StartValue=zeros(m,m);
maxd=size(A,2)/m-1;
options.NonZeroBlocks=0;

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% check whether G is known explicitly
G=MG1_EG(A,options.Verbose);
if (~isempty(G))
    return
end


numit=0;
check=1;
G=options.StartValue;
if (options.NonZeroBlocks == 0)
    % Shift Technique
    if (strfind(options.Mode,'Shift')>0)
        if (options.Verbose >= 1)
            Aold=A;
        end
        [A,drift,tau,v]=MG1_Shifts(A,options.ShiftType);
    end   
    if (strfind(options.Mode,'Natural')>0)
        while(check > 10^(-14) && numit < options.MaxNumIt)
            Gold=G;
            G=A(:,maxd*m+1:end);
            for j=maxd-1:-1:0
                G=A(:,j*m+1:(j+1)*m)+G*Gold;
            end
            check=norm(G-Gold,inf);
            numit=numit+1;
            if (~mod(numit,options.Verbose))
                fprintf('Check after %d iterations: %d\n',numit,check);
                drawnow;
            end
        end
    end

    if (strfind(options.Mode,'Traditional')>0)
        while(check > 10^(-14) && numit < options.MaxNumIt)
            Gold=G;
            G=A(:,maxd*m+1:end);
            for j=maxd-1:-1:2
                G=A(:,j*m+1:(j+1)*m)+G*Gold;
            end
            G=A(:,1:m)+G*Gold^2;
            G=(eye(m)-A(:,m+1:2*m))^(-1)*G;
            check=norm(G-Gold,inf);
            numit=numit+1;
            if (~mod(numit,options.Verbose))
                fprintf('Check after %d iterations: %d\n',numit,check);
                drawnow;
            end
        end
    end

    if (strfind(options.Mode,'U-Based')>0)
        while(check > 10^(-14) && numit < options.MaxNumIt)
            Gold=G;
            G=A(:,maxd*m+1:end);
            for j=maxd-1:-1:1
                G=A(:,j*m+1:(j+1)*m)+G*Gold;
            end
            G=(eye(m)-G)^(-1)*A(:,1:m);
            check=norm(G-Gold,inf);
            numit=numit+1;
            if (~mod(numit,options.Verbose))
                fprintf('Check after %d iterations: %d\n',numit,check);
                drawnow;
            end
        end
    end
    if (strfind(options.Mode,'Shift')>0)
        switch options.ShiftType
            case 'one'
                G=G+(drift<1)*ones(m,m)/m;
            case 'tau'
                G=G+(drift>1)*tau*v*ones(1,m);
            case 'dbl'
                G=G+(drift<1)*ones(m,m)/m+(drift>1)*tau*v*ones(1,m);
            end      
        if (options.Verbose >= 1)
            A=Aold;
        end    
    end
else
    if (strfind(options.Mode,'Shift')>0)
        warning('MATLAB:MG1_FI:ShiftIgnored',...
        'Shift is ignored due to NonZeroBlocks option');
    end    
    % Check length of vec option
    if (max(size(options.NonZeroBlocks) ~= [1 maxd]))
        error('MATLAB:MG1_FI:InvalidNonZeroBlockVec',...
                'The NonZeroBlocks option vector must be of length %d',maxd);
    end    
    vec=[0 options.NonZeroBlocks];
    vec=vec(2:end)-vec(1:end-1);
    % Check whether vec is increasing
    if (min(vec) < 1)
        error('MATLAB:MG1_FI:InvalidNonZeroBlockVec',...
                'The NonZeroBlocks option vector must be strictly increasing');
    end
    % If A1 = 0 then Traditional = Natural
    if (vec(1) > 1 && strcmp(options.Mode,'Traditional'))
        options.Mode='Natural';
    end
    if (strfind(options.Mode,'Natural')>0)
        while(check > 10^(-14) && numit < options.MaxNumIt)
            Gold=G;
            G=A(:,maxd*m+1:end);
            for j=maxd-1:-1:0
                G=A(:,j*m+1:(j+1)*m)+G*Gold^vec(j+1);
            end
            check=norm(G-Gold,inf);
            numit=numit+1;
            if (~mod(numit,options.Verbose))
                fprintf('Check after %d iterations: %d\n',numit,check);
                drawnow;
            end
        end
    end

    if (strfind(options.Mode,'Traditional')>0)
        % vec(1) = 1
        while(check > 10^(-14) && numit < options.MaxNumIt)
            Gold=G;
            G=A(:,maxd*m+1:end);
            for j=maxd-1:-1:2
                G=A(:,j*m+1:(j+1)*m)+G*Gold^vec(j+1);
            end
            G=A(:,1:m)+G*Gold^(1+vec(2));
            G=(eye(m)-A(:,m+1:2*m))^(-1)*G;
            check=norm(G-Gold,inf);
            numit=numit+1;
            if (~mod(numit,options.Verbose))
                fprintf('Check after %d iterations: %d\n',numit,check);
                drawnow;
            end
        end
    end

    if (strfind(options.Mode,'U-Based')>0)
        while(check > 10^(-14) && numit < options.MaxNumIt)
            Gold=G;
            G=A(:,maxd*m+1:end);
            for j=maxd-1:-1:1
                G=A(:,j*m+1:(j+1)*m)+G*Gold^vec(j+1);
            end
            G=(eye(m)-G*Gold^(vec(1)-1))^(-1)*A(:,1:m);
            check=norm(G-Gold,inf);
            numit=numit+1;
            if (~mod(numit,options.Verbose))
                fprintf('Check after %d iterations: %d\n',numit,check);
                drawnow;
            end
        end
    end
    
end
if (numit == options.MaxNumIt)
    warning('Maximum Number of Iterations %d reached',numit);
end

if (options.Verbose>0)
    Gcheck=A(:,maxd*m+1:end);
    for j=maxd-1:-1:0
        if (options.NonZeroBlocks == 0)
            Gcheck=A(:,j*m+1:(j+1)*m)+Gcheck*G;
        else
            Gcheck=A(:,j*m+1:(j+1)*m)+Gcheck*G^vec(j+1);
        end    
    end
    res_norm=norm(G-Gcheck,inf);
    fprintf('Final Residual Error for G: %d\n',res_norm);
end
