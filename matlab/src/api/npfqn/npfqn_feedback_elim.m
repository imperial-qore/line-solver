function result = npfqn_feedback_elim(P, rho, varargin)
% NPFQN_FEEDBACK_ELIM Near-immediate feedback elimination for RQNA.
%
% RESULT = NPFQN_FEEDBACK_ELIM(P, RHO) computes, for each station of an open
% queueing network with routing matrix P and traffic intensities RHO, the
% probability that a departing customer returns to that station WITHOUT passing
% through a busier one, and the modified service description that eliminates
% that feedback.
%
% WHY FEEDBACK BREAKS DECOMPOSITION. A parametric decomposition treats the
% arrival stream at each station as if it were renewal. Feedback destroys that
% badly: a customer that leaves a busy station and comes straight back arrives
% exactly when the station is busy, so the flow is strongly correlated with the
% queue it feeds. The fix is not to model the correlation but to REMOVE the
% feedback, by folding the repeated visits into the service time.
%
% THE TRANSFORMATION. With feedback probability p, a customer is served a
% geometric number of times, so the effective service is S_p = sum_{i<=N} S_i
% with N geometric of mean 1/(1-p). Hence
%   - effective mean service E[S]/(1-p),
%   - effective service SCV p + (1-p)cs^2 (eq. 37 and the line after it),
%   - fresh arrival rate lambda(1-p),
%   - per-visit waiting time = (1-p) times the wait in the modified system.
% The modified system has the SAME heavy-traffic limits for queue length,
% workload, waiting time and external departures, so the elimination is
% asymptotically exact rather than merely plausible.
%
% NEAR-IMMEDIATE, NOT JUST IMMEDIATE. Feedback rarely returns a customer in one
% hop. What matters is whether it returns WITHOUT PASSING A BUSIER STATION: a
% detour through a station of lower traffic intensity is fast on the time scale
% of the busy station, so it behaves like immediate feedback.
%
% Options:
%   'cs2', VEC          - service SCV per station; adds modifiedScv to the output
%   'lambda', VEC       - arrival rate per station; adds modifiedRates
%   'immediateOnly', TF - keep only the self-loops P(i,i), i.e. immediate
%                         feedback in the strict sense of Section 4.1
%
% Returns a struct with fields feedbackProb (p-hat per station), visitInflation
% (1/(1-p), the mean visits per customer), modifiedScv, modifiedRates,
% modifiedRouting and reductionExact.
%
% Example:
%   % the three-station example of Dai et al. (1994), Section 6.1 of the paper
%   P = [0 1 0; 0.5 0 0.5; 0 0.5 0];
%   res = npfqn_feedback_elim(P, [0.675 0.9 0.45]);
%   res.feedbackProb     % 0, 0.75, 0
%
% Reference: W. Whitt, W. You (2022). A robust queueing network analyzer based
% on indices of dispersion. Naval Research Logistics 69(1), 36-56, Section 4.
%
% See also NPFQN_TRAFFIC_IDC, SOLVER_RQNA.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

options = struct('cs2', [], 'lambda', [], 'immediateonly', false);
for i = 1:2:numel(varargin)
    if i+1 > numel(varargin)
        line_error(mfilename, sprintf('option %s has no value', char(varargin{i})));
    end
    name = lower(char(varargin{i}));
    if ~isfield(options, name)
        line_error(mfilename, sprintf('unknown option %s', char(varargin{i})));
    end
    options.(name) = varargin{i+1};
end

m = size(P, 1);
if size(P, 2) ~= m
    line_error(mfilename, 'The routing matrix must be square.');
end
if any(P(:) < -1e-12) || any(sum(P, 2) > 1 + 1e-9)
    line_error(mfilename, 'The routing matrix must be substochastic.');
end
rho = rho(:).';
if numel(rho) ~= m
    line_error(mfilename, 'One traffic intensity per station is required.');
end

phat = zeros(1, m);
for i = 1:m
    if options.immediateonly
        phat(i) = P(i,i);
        continue
    end
    % Stations a customer may pass through on a near-immediate return: those
    % NOT MORE loaded than i. A detour through a busier station is not fast on
    % the time scale of station i, so it is not near-immediate; one through a
    % station of equal load is, which is why the test is <= and not <. This is
    % the cloud of eqs. (3.8)-(3.9) with H = {i}, and the same one
    % solver_rqna.m applies -- the two must not drift.
    idx = setdiff(find(rho <= rho(i) + 1e-9), i);
    ret = P(i,i);
    if ~isempty(idx)
        % (I-Q)^-1 r is the probability of eventually reaching i from each
        % allowed station without leaving the allowed set.
        reach = (eye(numel(idx)) - P(idx,idx)) \ P(idx,i);
        ret = ret + P(i,idx)*reach;
    end
    phat(i) = min(max(ret, 0), 1 - 1e-12);
end

result.feedbackProb = phat;
result.visitInflation = 1 ./ (1 - phat);
if ~isempty(options.cs2)
    cs2 = options.cs2(:).';
    if numel(cs2) ~= m
        line_error(mfilename, 'One service SCV per station is required.');
    end
    result.modifiedScv = phat + (1 - phat).*cs2;      % eq. (37)
end
if ~isempty(options.lambda)
    lam = options.lambda(:).';
    if numel(lam) ~= m
        line_error(mfilename, 'One arrival rate per station is required.');
    end
    result.modifiedRates = lam .* (1 - phat);
end

% The reduced network. For IMMEDIATE feedback the reduction is exact and
% unambiguous: drop the self-loop and renormalize the rest of the row, since a
% customer that does not feed back goes where it would have gone anyway. For
% NEAR-IMMEDIATE feedback the return path runs through other stations, so there
% is no such row-local reduction; the elimination then applies to the SERVICE
% description at the station, which modifiedScv and modifiedRates carry.
Pmod = P;
for i = 1:m
    loop = P(i,i);
    if loop <= 0
        continue
    end
    Pmod(i,i) = 0;
    rest = sum(Pmod(i,:));
    if rest > 0
        Pmod(i,:) = Pmod(i,:) * (sum(P(i,:)) - loop)/rest;
    end
end
result.modifiedRouting = Pmod;
result.reductionExact = options.immediateonly || all(abs(phat - diag(P)') < 1e-12);
end
