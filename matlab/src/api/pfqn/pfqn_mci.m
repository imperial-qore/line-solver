%{
%{
 % @file pfqn_mci.m
 % @brief Monte Carlo Integration (MCI) for normalizing constant.
%}
%}

%{
%{
 % @brief Monte Carlo Integration (MCI) for normalizing constant.
 % @fn pfqn_mci(D, N, Z, I, variant)
 % @param D Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param I Number of samples (default: 1e5).
 % @param variant MCI variant ('mci', 'imci', 'rm', default: 'imci').
 % @return G Normalizing constant estimate.
 % @return lG Logarithm of normalizing constant.
 % @return lZ Individual random sample log values.
%}
%}
function [G,lG,lZ] = pfqn_mci(D,N,Z,I,variant)
% [G,LG,LZ] = PFQN_MCI(D,N,Z,I,VARIANT)
%
% Normalizing constant estimation via Monte Carlo Integration
%
% Syntax:
% [G,lG,lZ] = pfqn_mci(D,N,Z,I,VARIANT)
% Input:
% D - demands (queues x classes)
% N - populations (1 x classes)
% Z - think times (1 x classes)
% I - samples
% VARIANT - 'mci', 'imci', 'rm'
%
% Output:
% lG - estimate of logG
% lZ - individual random samples
%
% Note: if the script returns a floating point range exception,
% double(log(mean(exp(sym(lZ))))) provides a better estimate of lG, but it
% is very time consuming due to the symbolic operations.
%
% Implementation: Giuliano Casale (g.casale@imperial.ac.uk), 16-Aug-2013

if nargin<3%~exist('Z','var')
    Z=0*N;
end
if nargin<4%~exist('I','var')
    I=1e5;
end
if nargin<5%~exist('variant','var')
    variant='imci';
end

[M,R] = size(D);

if isempty(D) || sum(D(:))<1e-4
    lGn = - sum(factln(N)) + sum(N.*log(sum(Z,1)));
    G=exp(lGn);
    lZ=[];
    return
end

%tput = N./(Z+sum(D)+max(D).*(sum(N)-1)); % balanced job bounds
%tput = N./Z; % balanced job bounds

%% IMCI
if strcmpi(variant,'imci') % improved mci
    tput = pfqn_bs(D,N,Z);
    util = D*tput';
    gamma = max( 0.01, 1-util )'; % MonteQueue 2.0 recommendation
elseif strcmpi(variant,'mci') % original mci
    tput = pfqn_bs(D,N,Z);
    util = D*tput';
    %% Original MCI
    for i=1:length(util)
        if util(i)>0.9
            gamma(i) = 1/sqrt(max(N)); % MonteQueue 2.0 recommendation
        else
            gamma(i) = 1-util(i); % MonteQueue 2.0 recommendation
        end
    end
elseif strcmpi(variant,'rm') % repairman problem
    tput = N./(sum(D,1)+Z+max(D,1)*(sum(N)-1)); % a single queue
    util = D*tput';
    %% Original MCI
    for i=1:length(util)
        if util(i)>0.9
            gamma(i) = 1/sqrt(max(N)); % MonteQueue 2.0 recommendation
        else
            gamma(i) = 1-util(i); % MonteQueue 2.0 recommendation
        end
    end
end
try
    for r=1:R
        logfact(r) = sum(log(1:N(r)));  % log N(r)!
    end
    
    % uniform sampling
    persistent VL
    if isempty(VL) || size(VL,1) < I || size(VL,2) < M
        VL = log(rand(I,M));
    end
    V = repmat(-1./gamma,I,1).*VL(1:I,1:M);
    ZI = repmat(Z,I,1);
    % importance sampling
    lZ = -(ones(1,M) - gamma) * V' - sum(log(gamma)) - sum(logfact) + N*log(V*D+ZI)';
    
    
    lG = logmeanexp(lZ); % return average    
    if isinf(lG)
        %    line_warning(mfilename,'Floating-point range exception, Monte Carlo integration will return an approximation.');
        lG = max(lZ);
    end
    G=exp(lG);    
catch ME
    getReport(ME,'basic')
end
end
