function G=MG1_CR(D,varargin)
%MG1_CR Cyclic reduction for M/G/1-Type Markov Chains [Bini,Meini] 
%
%   G=MG1_CR(A) computes the minimal nonnegative solution to the 
%   matrix equation G = A0 + A1 G + A2 G^2 + A3 G^3 + ... + A_max G^max, 
%   where A = [A0 A1 A2 A3 ... A_max] has m rows and m*max columns and is
%   a nonnegative matrix, with (A0+A1+A2+...+A_max) irreducible and 
%   stochastic
%
%   Optional Parameters:
%   
%       MaxNumIt: Maximum number of iterations (default: 50)
%       MaxNumRoot: Maximum number of roots used by Point-Wise CR (default: 2048)
%       EpsilonValue: Required accuracy at each step (default: 10^(-16)) 
%       Verbose: The residual error is printed at each step when set to 1,  
%                (default:0)
%       ShiftType: 'one'
%                  'tau'
%                  'dbl'
%       Mode: 'PWCR' uses the Point-Wise Cyclic Reduction  
%             'ShiftPWCR' uses the Shift + Point-Wise Cyclic Reduction 
%             (default: 'ShiftPWCR' with ShiftType='one')

OptionNames=['Mode        '; 
             'MaxNumIt    ';
             'ShiftType   ';
             'MaxNumRoot  ';
             'EpsilonValue';
             'Verbose     '];
OptionTypes=['char   '; 
             'numeric';
             'char   ';
             'numeric';
             'numeric';
             'numeric'];
OptionValues{1}=['ShiftPWCR';
                 'PWCR     '];

OptionValues{3}=['one';
                 'tau';
                 'dbl'];

             
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
%options.ProgressBar=0;
options.Mode='ShiftPWCR';
options.ShiftType='one';
options.MaxNumRoot=2048;
options.MaxNumIt=50;
options.Verbose=0;
options.EpsilonValue=10^(-16);

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% check whether G is known explicitly
G=MG1_EG(D,options.Verbose);
if (~isempty(G))
    return
end

if (strfind(options.Mode,'ShiftPWCR')>0)
    if (options.Verbose == 1)
        Dold=D;
    end
    [D,drift,tau,v]=MG1_Shifts(D,options.ShiftType);  
end

% start Cyclic Reduction
lastiter=0;
m=size(D,1);
D=D';
D=[D; zeros((2^(1+floor(log2(size(D,1)/m-1)))+1)*m-size(D,1),m)];

% Step 0
G=zeros(m,m);
Aeven=D(find(mod(kron(1:size(D,1)/m,ones(1,m)),2)),:);
Aodd=D(find(~mod(kron(1:size(D,1)/m,ones(1,m)),2)),:);

Ahatodd=[Aeven(m+1:end,:); D(end-m+1:end,:)];
Ahateven=Aodd;

Rj=D(m+1:2*m,:);
for i=3:size(D,1)/m
    Rj=Rj+D((i-1)*m+1:i*m,:);
end
Rj=inv(eye(m)-Rj);
Rj=D(1:m,:)*Rj;

numit=0;

while(1 && numit < options.MaxNumIt)
    numit=numit+1;
    nj=size(Aodd,1)/m-1;
    if (nj > 0)
        % Evaluate the 4 functions in the nj+1 roots using FFT
        % prepare for FFTs (such that they can be performed in 4 calls)
        temp1=reshape(Aodd(1:(nj+1)*m,:)',m^2,(nj+1))';
        temp2=reshape(Aeven(1:(nj+1)*m,:)',m^2,(nj+1))';
        temp3=reshape(Ahatodd(1:(nj+1)*m,:)',m^2,(nj+1))';
        temp4=reshape(Ahateven(1:(nj+1)*m,:)',m^2,(nj+1))';

        % FFTs
        temp1=fft(temp1,(nj+1));
        temp2=fft(temp2,(nj+1));
        temp3=fft(temp3,(nj+1));
        temp4=fft(temp4,(nj+1));

        % reform the 4*(nj+1) matrices
        temp1=reshape(temp1.',m,m*(nj+1)).';
        temp2=reshape(temp2.',m,m*(nj+1)).';
        temp3=reshape(temp3.',m,m*(nj+1)).';
        temp4=reshape(temp4.',m,m*(nj+1)).';

        % Next, we perform a point-wise evaluation of (6.20) - Thesis Meini
        for cnt=1:nj+1
            Ahatnew((cnt-1)*m+1:cnt*m,1:m)=temp4((cnt-1)*m+1:cnt*m,:)+...
                temp2((cnt-1)*m+1:cnt*m,:)*inv(eye(m)-temp1((cnt-1)*m+1:cnt*m,:))*temp3((cnt-1)*m+1:cnt*m,:);
            Anew((cnt-1)*m+1:cnt*m,1:m)=exp(-(cnt-1)*2j*pi/(nj+1))*temp1((cnt-1)*m+1:cnt*m,:)+...
                temp2((cnt-1)*m+1:cnt*m,:)*inv(eye(m)-temp1((cnt-1)*m+1:cnt*m,:))*temp2((cnt-1)*m+1:cnt*m,:);
        end

        % We now invert the FFTs to get Pz and Phatz

        % prepare for IFFTs (in 2 calls)
        Ahatnew=reshape(Ahatnew(1:(nj+1)*m,:).',m^2,(nj+1)).';
        Anew=reshape(Anew(1:(nj+1)*m,:).',m^2,(nj+1)).';

        % IFFTs
        Ahatnew=real(ifft(Ahatnew,(nj+1)));
        Anew=real(ifft(Anew,(nj+1)));

        % reform matrices Pi and Phati
        Ahatnew=reshape(Ahatnew',m,m*(nj+1))';
        Anew=reshape(Anew',m,m*(nj+1))';
    else % series Aeven, Aodd, Ahateven and Ahatodd are constant
        temp=Aeven*(eye(m)-Aodd)^(-1);
        Ahatnew=Ahateven+temp*Ahatodd;
        Anew=[temp*Aeven; Aodd];
        Aodd1=Aodd;
    end

    nAnew=0;
    deg=size(Anew,1)/m;
    for i=deg/2:deg-1
        nAnew=max(nAnew,norm(Anew(i*m+1:(i+1)*m,:),inf));
    end
    nAhatnew=0;
    deghat=size(Ahatnew,1)/m;
    for i=deghat/2:deghat-1
        nAhatnew=max(nAhatnew,norm(Ahatnew(i*m+1:(i+1)*m,:),inf));
    end 
    
    % c) the test
    while ((nAnew > (nj+1)*options.EpsilonValue || ...
            (nAhatnew > (nj+1)*options.EpsilonValue)) && ...
            nj+1 < options.MaxNumRoot)

        nj=2*(nj+1)-1;
        stopv=min([nj+1 size(Aodd,1)/m]);

        % prepare for FFTs
        temp1=reshape(Aodd(1:stopv*m,:)',m^2,stopv)';
        temp2=reshape(Aeven(1:stopv*m,:)',m^2,stopv)';
        temp3=reshape(Ahatodd(1:stopv*m,:)',m^2,stopv)';
        temp4=reshape(Ahateven(1:stopv*m,:)',m^2,stopv)';

        % FFTs
        temp1=fft(temp1,(nj+1),1);
        temp2=fft(temp2,(nj+1),1);
        temp3=fft(temp3,(nj+1),1);
        temp4=fft(temp4,(nj+1),1);

        % reform the 4*(nj+1) matrices
        temp1=reshape(temp1.',m,m*(nj+1)).';
        temp2=reshape(temp2.',m,m*(nj+1)).';
        temp3=reshape(temp3.',m,m*(nj+1)).';
        temp4=reshape(temp4.',m,m*(nj+1)).';

        % Next, we perform a point-wise evaluation of (6.20) - Thesis Meini
        for cnt=1:nj+1
            Ahatnew((cnt-1)*m+1:cnt*m,1:m)=temp4((cnt-1)*m+1:cnt*m,:)+...
                temp2((cnt-1)*m+1:cnt*m,:)*inv(eye(m)-temp1((cnt-1)*m+1:cnt*m,:))*temp3((cnt-1)*m+1:cnt*m,:);
            Anew((cnt-1)*m+1:cnt*m,1:m)=exp(-(cnt-1)*2j*pi/(nj+1))*temp1((cnt-1)*m+1:cnt*m,:)+...
                temp2((cnt-1)*m+1:cnt*m,:)*inv(eye(m)-temp1((cnt-1)*m+1:cnt*m,:))*temp2((cnt-1)*m+1:cnt*m,:);
        end

        % We now invert the FFTs to get Pz and Phatz
        % prepare for IFFTs
        Ahatnew=reshape(Ahatnew(1:(nj+1)*m,:).',m^2,(nj+1)).';
        Anew=reshape(Anew(1:(nj+1)*m,:).',m^2,(nj+1)).';

        % IFFTs
        Ahatnew=real(ifft(Ahatnew,(nj+1)));
        Anew=real(ifft(Anew,(nj+1)));

        % reform matrices Pi and Phati
        Ahatnew=reshape(Ahatnew',m,m*(nj+1))';
        Anew=reshape(Anew',m,m*(nj+1))';

        vec1=zeros(1,m);
        vec2=zeros(1,m);
        for i=1:size(Anew,1)/m-1
            vec1=vec1+i*sum(Anew(i*m+1:(i+1)*m,:));
            vec2=vec2+i*sum(Ahatnew(i*m+1:(i+1)*m,:));
        end
        nAnew=0;
        deg=size(Anew,1)/m;
        for i=deg/2:deg-1
            nAnew=max(nAnew,norm(Anew(i*m+1:(i+1)*m,:),inf));
        end
        nAhatnew=0;
        deghat=size(Ahatnew,1)/m;
        for i=deghat/2:deghat-1
            nAhatnew=max(nAhatnew,norm(Ahatnew(i*m+1:(i+1)*m,:),inf));
        end
    end
    if (( nAnew > (nj+1)*options.EpsilonValue || ...
            nAhatnew > (nj+1)*options.EpsilonValue ) && ...
            nj+1 >= options.MaxNumRoot)
        warning('MATLAB:MG1_CR:MaxNumRootReached',...
            'Maximum number of ''%d'' reached, accuracy might be affected',options.MaxNumRoot);
    end

    if (nj > 1)
        Anew=Anew(1:m*(nj+1)/2,:);
        Ahatnew=Ahatnew(1:m*(nj+1)/2,:);
    end
    
    % compute Aodd, Aeven, ...
    Aeven=Anew(find(mod(kron(1:size(Anew,1)/m,ones(1,m)),2)),:);
    Aodd=Anew(find(~mod(kron(1:size(Anew,1)/m,ones(1,m)),2)),:);
    
    Ahateven=Ahatnew(find(mod(kron(1:size(Ahatnew,1)/m,ones(1,m)),2)),:);
    Ahatodd=Ahatnew(find(~mod(kron(1:size(Ahatnew,1)/m,ones(1,m)),2)),:);
    
    if (options.Verbose==1)
        if (strcmp(options.Mode,'PWCR'))
            fprintf('The Point-wise evaluation of Iteration %d required %d roots\n',numit,nj+1);
        else
            fprintf('The Shifted PWCR evaluation of Iteration %d required %d roots\n',numit,nj+1);
        end
        drawnow;
    end
    
    % test stopcriteria
    if (strcmp(options.Mode,'PWCR') || strcmp(options.Mode,'DCR'))
        Rnewj=Anew(m+1:2*m,:);
        for i=3:size(Anew,1)/m
            Rnewj=Rnewj+Anew((i-1)*m+1:i*m,:);
        end
        Rnewj=inv(eye(m)-Rnewj);
        Rnewj=Anew(1:m,:)*Rnewj;
        if (max(max(abs(Rj-Rnewj))) < options.EpsilonValue || ...
                max(sum(eye(m)-Anew(1:m,:)*inv(eye(m)-Anew(m+1:2*m,:)))) < options.EpsilonValue)
            G=Ahatnew(1:m,:);
            for i=2:size(Ahatnew,1)/m
                G=G+Rnewj*Ahatnew((i-1)*m+1:i*m,:);
            end
            G=D(1:m,:)*inv(eye(m)-G);
            break
        end
        Rj=Rnewj;
        % second condition tests whether Ahatnew is degree 0 (numerically)
        if (norm(Anew(1:m,1:m))<options.EpsilonValue || sum(sum(Ahatnew(m+1:end,:)))<options.EpsilonValue || max(sum(eye(m)-D(1:m,:)*inv(eye(m)-Ahatnew(1:m,:)))) < options.EpsilonValue)
            G=D(1:m,:)*inv(eye(m)-Ahatnew(1:m,:));
            break
        end
    else 
        Gold=G;
        G=D(1:m,:)*inv(eye(m)-Ahatnew(1:m,:));
        if (norm(G-Gold,inf)<options.EpsilonValue || norm(Ahatnew(m+1:end,:),inf)<options.EpsilonValue)
            break;
        end     
    end
end
if (numit == options.MaxNumIt & isempty(G))
    warning('Maximum Number of Iterations %d reached',numit);
    G=D(1:m,:)*inv(eye(m)-Ahatnew(1:m,:));
end    

G=G';

if (strcmp(options.Mode,'ShiftPWCR'))
    switch options.ShiftType
        case 'one'
            G=G+(drift<1)*ones(m,m)/m;
        case 'tau'
            G=G+(drift>1)*tau*v*ones(1,m);
        case 'dbl'
            G=G+(drift<1)*ones(m,m)/m+(drift>1)*tau*v*ones(1,m);
    end    
end    
    
if (options.Verbose==1)
    if (strcmp(options.Mode,'PWCR'))
        D=D';
    else
        D=Dold;
    end
    temp=D(:,end-m+1:end);
    for i=size(D,2)/m-1:-1:1
        temp=D(:,(i-1)*m+1:i*m)+temp*G;
    end
    res_norm=norm(G-temp,inf);
    fprintf('Final Residual Error for G: %d\n',res_norm);
end
