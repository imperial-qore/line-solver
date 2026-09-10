function tf = mam_is_renewal_map(D0, D1)
% MAM_IS_RENEWAL_MAP  True if the process (D0,D1) is a renewal process.
%
% A renewal process embedded as a MAP has D1 = t * alpha (rank one): the phase
% entered after an event does not depend on the phase the process was in at that
% event, so successive interarrival times are independent. Equivalently, the
% process carries no autocorrelation and is fully described by its interevent
% marginal.
%
% The test is on the algebraic structure alone, so it holds equally for a MAP,
% for a matrix-exponential (ME) and for a rational arrival process (RAP): an ME
% renewal stream satisfies it, a correlated MAP or RAP does not. That is what
% makes it the right guard in front of any closed form that reads only the
% marginal (a PH pair (pie, D0), a mean and an SCV), because such a form is
% exact for a renewal input and silently discards the correlation otherwise.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

ns = size(D0, 1);
if ns == 1
    tf = true;   % a one-phase process is memoryless, hence renewal
    return;
end
tExit = D1 * ones(ns, 1);
alpha = map_pie({D0, D1});
tf = norm(D1 - tExit * alpha, 'fro') < 1e-9 * max(1, norm(D1, 'fro'));
end
