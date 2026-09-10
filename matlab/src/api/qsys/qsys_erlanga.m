function result = qsys_erlanga(lambda, mu, theta, s, r, varargin)
% QSYS_ERLANGA Exact analysis of the Erlang A model M/M/s/r+M.
%
% RESULT = QSYS_ERLANGA(LAMBDA, MU, THETA, S) analyzes the Erlang A queue:
% Poisson arrivals at rate LAMBDA, exponential service of rate MU at each of S
% servers, and exponential patience of rate THETA, so a waiting customer
% abandons after an exponential time of mean 1/THETA. The waiting room is
% infinite; THETA > 0 makes the model ergodic at every load, including LAMBDA
% above S*MU.
%
% RESULT = QSYS_ERLANGA(LAMBDA, MU, THETA, S, R) allows only R extra waiting
% spaces, so an arrival finding S+R customers is blocked and lost. R = Inf is
% the default.
%
% The number in system is the birth-and-death process with birth rate LAMBDA and
% death rate min(k,S)*MU + (k-S)^+ *THETA, so every measure below is EXACT: this
% is the special case in which the state-dependent Markovian approximation of
% QSYS_MGISRGI_WHITT reproduces the model rather than approximating it (eq. 7.12
% of the reference). THETA = 0 recovers M/M/s/r without abandonment, and then a
% finite R is required whenever LAMBDA >= S*MU.
%
% Accepts and returns exactly what QSYS_MGISRGI_WHITT does, including the
% 'wPoints' option for the waiting-time cdfs.
%
% Example:
%   res = qsys_erlanga(102, 1/10, 1/10, 100);
%   res.probAbandon    % fraction of arrivals that abandon
%
% Reference: W. Whitt (2005). Engineering solution of a basic call-center model.
% Management Science 51(2), 221-235, Section 7 and eq. (7.12). The model itself
% is due to C. Palm (1937, 1957).
%
% See also QSYS_MGISRGI_WHITT, QSYS_MMK, QSYS_MMCK.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 5 || isempty(r)
    r = Inf;
end
if theta <= 0 && ~isfinite(r) && lambda >= s*mu
    line_error(mfilename, ['without abandonment (theta = 0) and with an infinite waiting ' ...
        'room the queue is unstable at lambda >= s*mu; give a finite r or a positive theta']);
end

result = qsys_mgisrgi_whitt(lambda, mu, s, r, theta, varargin{:});
end
