function varargout=qsys_mginf(lambda,mu,varargin)
% [L,Lq,W,Wq,p0]=QSYS_MGINF(LAMBDA,MU)
% [L,Lq,W,Wq,p0,pk]=QSYS_MGINF(LAMBDA,MU,K)
%
% Exact solution for M/G/∞ queue (infinite servers).
% Performance is independent of the service time distribution shape (G).
% The number of customers in the system follows a Poisson distribution.
%
% Input:
%   lambda - arrival rate
%   mu - service rate (mean service time = 1/mu)
%   K (optional) - state for probability computation
%
% Output:
%   L - average number of customers in system
%   Lq - average number of customers in queue (always 0)
%   W - average time in system (= 1/mu)
%   Wq - average waiting time in queue (always 0)
%   p0 - probability of empty system
%   pk (optional) - probability of exactly k customers in system

rho = lambda / mu;

% Exact results for M/G/∞
L = rho;           % Average number of customers = average busy servers
Lq = 0;            % No queue since infinite servers
W = 1 / mu;        % Average time = service time (independent of arrival process)
Wq = 0;            % No waiting time
p0 = exp(-rho);    % Empty system probability (Poisson distribution)

% Always return 5 outputs
varargout{1} = L;
varargout{2} = Lq;
varargout{3} = W;
varargout{4} = Wq;
varargout{5} = p0;

% If k is provided, compute and return P(k customers in system)
if nargin >= 3
    k = varargin{1};
    pk = exp(-rho) * rho^k / factorial(k);
    varargout{6} = pk;
end
end
