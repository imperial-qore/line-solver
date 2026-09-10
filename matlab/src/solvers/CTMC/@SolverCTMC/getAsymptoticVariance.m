function result = getAsymptoticVariance(self, f)
% RESULT = GETASYMPTOTICVARIANCE(F)
%
% The asymptotic variance of the time-average of a reward F along a sample path
% of this model's CTMC.
%
% WHAT IT IS FOR. A simulation estimate of a steady-state mean has a standard
% error that shrinks like sqrt(sigma^2/t), where sigma^2 is NOT the stationary
% variance of F but its ASYMPTOTIC variance, which also carries the
% autocorrelation of the path. That number is what says how long a run has to
% be, and SIM_RUNLENGTH turns it into a run length for a target precision. It
% cannot be guessed from the stationary variance: on M/M/1 the two differ by a
% factor that blows up like (1-rho)^-2.
%
% F is either a vector of one reward value per CTMC state, in the state order
% GETGENERATOR returns, or a function handle applied to each row of the state
% space.
%
% Returns the struct of SIM_ASYMVAR_CTMC: mean, variance, asymptoticVariance and
% the deviation vector d that solves A d = -(F - E_pi[F]).
%
% Example:
%   sol = SolverCTMC(model);
%   res = sol.getAsymptoticVariance(@(state) state(2));   % queue length
%   plan = sim_runlength(res.mean, res.asymptoticVariance, 'relprecision', 0.05);
%
% Reference: W. Whitt (1989). Planning queueing simulations. Management Science
% 35(11), 1341-1366.
%
% See also SIM_ASYMVAR_CTMC, SIM_RUNLENGTH, GETGENERATOR.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

infGen = self.getGenerator();
if issparse(infGen)
    infGen = full(infGen);
end
n = size(infGen, 1);

if isa(f, 'function_handle')
    % One reward per state, evaluated on the state space the generator was
    % built from: SolverCTMC.getStateSpace returns it in the same order.
    space = self.getStateSpace();
    if size(space, 1) ~= n
        line_error(mfilename, sprintf(['the state space has %d rows but the generator is %dx%d; ' ...
            'pass the reward as a vector instead'], size(space,1), n, n));
    end
    fvec = zeros(n, 1);
    for i = 1:n
        fvec(i) = f(space(i,:));
    end
else
    fvec = f(:);
    if numel(fvec) ~= n
        line_error(mfilename, sprintf(['the reward vector has %d entries but the generator is ' ...
            '%dx%d'], numel(fvec), n, n));
    end
end

pi_ss = ctmc_solve(infGen);
result = sim_asymvar_ctmc(infGen, fvec, pi_ss(:));
end
