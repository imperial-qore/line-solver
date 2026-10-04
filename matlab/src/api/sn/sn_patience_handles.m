function h = sn_patience_handles(sn, ist, r)
% H = SN_PATIENCE_HANDLES(SN, IST, R)
%
% Patience (time-to-abandon) handles for station IST, class R.
%
% The abandonment solvers need the patience law as FUNCTIONS -- a complementary
% cdf and a hazard rate -- not as moments, because that is what the underlying
% theory consumes: Whitt's engineering solution reads the hazard near the
% origin, and the fluid models integrate the ccdf. LINE stores the law as a
% MAP/PH pair in sn.impatienceProc, from which both are available in closed
% form.
%
% Returns [] when the station-class pair has no reneging patience configured,
% otherwise a struct with fields ccdf, pdf, hazard (function handles), mean,
% isExponential and rate.
%
% See also: python/line_solver/api/sn/patience.py

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

h = [];
if ~isfield(sn, 'impatienceClass') || isempty(sn.impatienceClass)
    return
end
if size(sn.impatienceClass,1) < ist || size(sn.impatienceClass,2) < r
    return
end
if sn.impatienceClass(ist, r) ~= ImpatienceType.RENEGING
    return
end

rate = 0;
if isfield(sn, 'impatienceMu') && ~isempty(sn.impatienceMu)
    rate = sn.impatienceMu(ist, r);
end
isExp = false;
if isfield(sn, 'impatienceType') && ~isempty(sn.impatienceType)
    isExp = sn.impatienceType(ist, r) == ProcessType.EXP;
end

pair = [];
if isfield(sn, 'impatienceProc') && ~isempty(sn.impatienceProc) ...
        && size(sn.impatienceProc,1) >= ist && size(sn.impatienceProc,2) >= r
    pair = sn.impatienceProc{ist, r};
end

if isempty(pair)
    if rate <= 0
        return
    end
    % Only the rate is on record, so the law is exponential by construction.
    h = struct('ccdf', @(t) exp(-rate*t), ...
        'pdf', @(t) rate*exp(-rate*t), ...
        'hazard', @(t) rate*ones(size(t)), ...
        'mean', 1/rate, ...
        'isExponential', true, ...
        'rate', rate);
    return
end

D0 = pair{1};
D1 = pair{2};
ccdf = @(t) 1 - map_cdf({D0,D1}, t);
pdf = @(t) map_pdf({D0,D1}, t);
if rate > 0
    meanval = 1/rate;
else
    meanval = Inf;
end
h = struct('ccdf', ccdf, ...
    'pdf', pdf, ...
    'hazard', @(t) local_hazard(ccdf, pdf, t, rate), ...
    'mean', meanval, ...
    'isExponential', isExp, ...
    'rate', rate);
end

function v = local_hazard(ccdf, pdf, t, rate)
% h = f/(1-F). Past the point where the ccdf underflows the hazard is the
% asymptotic decay rate, and returning that is better conditioned than
% dividing two zeros.
c = ccdf(t);
v = zeros(size(c));
ok = c > 1e-300;
f = pdf(t);
v(ok) = f(ok)./c(ok);
if any(~ok(:))
    if rate > 0
        v(~ok) = rate;
    else
        v(~ok) = 0;
    end
end
end
