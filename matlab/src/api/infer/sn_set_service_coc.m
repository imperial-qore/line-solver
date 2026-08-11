function sn = sn_set_service_coc(sn, stationIdx, classIdx, rate, scv)
% SN_SET_SERVICE_COC Update service rate preserving cell-of-cells format.
%
% Works like sn_set_service with autoRefresh=true, but writes mu, phi,
% proc, pie using cell-of-cells indexing (sn.mu{i}{k}) instead of 2D
% cell indexing (sn.mu{i,k}). This avoids reshaping the cell arrays,
% which would break solvers that expect cell-of-cells format.
%
% Inputs:
%   sn         - NetworkStruct
%   stationIdx - station index (1-based)
%   classIdx   - class index (1-based)
%   rate       - new service rate (positive scalar)
%   scv        - squared coefficient of variation (default 1.0)
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
% This code is released under the 3-Clause BSD License.

if nargin < 5 || isempty(scv)
    scv = 1.0;
end

% Update rates and scv matrices
sn.rates(stationIdx, classIdx) = rate;
sn.scv(stationIdx, classIdx) = scv;

% Skip if rate is invalid
if isnan(rate) || rate <= 0 || isinf(rate)
    return;
end

meanVal = 1.0 / rate;

% Determine process type and create MAP based on SCV
if isnan(scv) || abs(scv - 1.0) < 1e-10
    MAP = map_exponential(meanVal);
    nPhases = 1;
    procType = ProcessType.EXP;
elseif scv < 1.0
    k = max(1, ceil(1.0 / scv));
    MAP = map_erlang(meanVal, k);
    nPhases = k;
    procType = ProcessType.ERLANG;
else
    MAP = map_hyperexp(meanVal, scv);
    if ~isempty(MAP)
        nPhases = 2;
        procType = ProcessType.HYPEREXP;
    else
        MAP = map_exponential(meanVal);
        nPhases = 1;
        procType = ProcessType.EXP;
    end
end

D0 = MAP{1};
D1 = MAP{2};

% Update in cell-of-cells format
sn.proc{stationIdx}{classIdx} = MAP;
sn.procid(stationIdx, classIdx) = procType;
sn.phases(stationIdx, classIdx) = nPhases;
sn.phasessz(stationIdx, classIdx) = max(nPhases, 1);

% Recompute phaseshift for this station
cumSum = 0;
sn.phaseshift(stationIdx, 1) = 0;
for c = 1:sn.nclasses
    cumSum = cumSum + sn.phasessz(stationIdx, c);
    if c + 1 <= size(sn.phaseshift, 2)
        sn.phaseshift(stationIdx, c + 1) = cumSum;
    end
end

% Update mu (rates from -diag(D0))
sn.mu{stationIdx}{classIdx} = -diag(D0);

% Update phi (completion probabilities)
phiVec = zeros(nPhases, 1);
for i = 1:nPhases
    d1RowSum = sum(D1(i, :));
    d0Diag = -D0(i, i);
    if d0Diag ~= 0
        phiVec(i) = d1RowSum / d0Diag;
    end
end
sn.phi{stationIdx}{classIdx} = phiVec;

% Update pie (initial phase distribution)
sn.pie{stationIdx}{classIdx} = map_pie(MAP);

end
