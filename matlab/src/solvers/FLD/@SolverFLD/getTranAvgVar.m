function [t, QVart, Sigmat, QCovt] = getTranAvgVar(self)
% [T, QVART, SIGMAT, QCOVT] = GETTRANAVGVAR()
%
% Transient queue-length VARIANCE per station and class.
%
% Two methods compute a second moment along the trajectory. 'kp' integrates the
% covariance of the Ko-Pender diffusion limit alongside the fluid mean (see
% solver_fluid_kp.m). 'dae' integrates the linear-noise covariance alongside the
% min-normal mean as one differential-algebraic system (see solver_fluid_dae.m).
% Every other fluid method carries a first moment only -- 'minnormal' included,
% whose covariance is the STATIONARY one and does not vary along the horizon --
% so this errors rather than returning a misleading zero.
%
% QVART is an (nstations x nclasses) cell of column vectors over the time grid
% T. SIGMAT is the full state covariance, a (dim x dim x numel(T)) array, so
% cross-station and cross-class terms survive rather than only the per-block
% totals in QVART. QCOVT is SIGMAT aggregated onto station-class pairs, an
% (M*K) x (M*K) x numel(T) array indexed IR = (R-1)*M + I: its diagonal is
% QVART and its off-diagonal entries are the cross terms, in an index space a
% caller can use without knowing the method's own phase layout. SolverENV's
% 'meancov' coupling reads it and seeds the next stage through
% options.config.init_qlen / init_qcov, which are in the same index space.
%
% The lang='cpp' route does not carry QCOVT: line-cli reports the phase-layout
% SIGMAT only, and this function does not hold the block layout needed to
% aggregate it. It comes back empty there rather than as a misleading zero.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.


isdae = any(strcmp(self.options.method, {'dae','fld.dae'}));
if ~isdae && ~any(strcmp(self.options.method, {'kp','fld.kp'}))
    line_error(mfilename, ['getTranAvgVar needs options.method=''kp'' or ''dae''; the other fluid ' ...
        'methods integrate the mean only and carry no second moment.']);
end

options = self.getOptions;
if ~isdae
    options.method = 'kp';
end
% THE INITIAL CONDITION IS RESOLVED HERE, as SOLVER_FLUID_ANALYZER does before
% every other fluid route. Both functions below take their starting point from
% options.init_sol, and SOLVER_FLUID falls back to an even spread of the
% population over the stations when it is empty -- NOT to the model's own
% initial state. This function calls them directly rather than through the
% analyzer, so leaving it empty started every transient from that spread: the
% mean reported by getTranAvg honoured initFromMarginal and the second moment
% reported here did not, silently, and the two described different trajectories.
% SolverENV's 'meancov' coupling is what made it visible, since it hands each
% stage an entry state and reads both back.
if isempty(options.init_sol)
    options.init_sol = solver_fluid_initsol(self.model.getStruct(), options);
end
% Steady state and transient are mutually exclusive at the API level, so a
% caller that left the horizon unbounded gets one resolved the same way
% getTranAvg does, from the slowest rate in the model.
if ~isfinite(options.timespan(2))
    sn = self.model.getStruct();
    rates = sn.rates(isfinite(sn.rates) & sn.rates > 0);
    if isempty(rates)
        slow = 1;
    else
        slow = min(rates);
    end
    options.timespan = [0, max(10, 30/slow)];
end

% lang='cpp' takes the second moments from line-cli (-s fluid -a tranvar),
% which integrates the same kp fluid AND diffusion limits over the horizon
% resolved just above. It is handed the RESOLVED options, so both paths use
% the same slowest-rate rule when the caller left the timespan unbounded.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    [t, QVart, Sigmat] = CPPLINE.tranAvgVar(self.name, self.model, options);
    QCovt = [];
    return
end

if isdae
    % The DAE route reports its transient second moment through the moment
    % struct, on the SAME time grid ode15s produced for the mean, so the two
    % need no interpolation onto one another. An empty Sigmat means the
    % covariance exceeded options.config.dae_maxcov and was held stationary;
    % say so rather than returning zeros that look like a computed answer.
    [~,~,~,~,~,~,~,~,~,~,~,~,moments] = solver_fluid_dae(self.model.getStruct(), options);
    if isempty(moments.Sigmat)
        line_error(mfilename, ['The dae method held the covariance at its stationary value because ' ...
            'the model exceeds options.config.dae_maxcov, so there is no transient second moment ' ...
            'to report. Raise that limit, or read the stationary covariance with getMoments.']);
    end
    t = moments.tvar;
    QVart = moments.QVart;
    Sigmat = moments.Sigmat;
    QCovt = moments.QCovt;
    return
end

[~,~,~,~,~,~,~,~,~,t,~,~,QVart,Sigmat,QCovt] = solver_fluid_kp(self.model.getStruct(), options);
end
