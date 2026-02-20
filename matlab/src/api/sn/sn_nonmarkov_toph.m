function sn = sn_nonmarkov_toph(sn, options)
% SN = SN_NONMARKOV_TOPH(SN, OPTIONS)
% Convert non-Markovian distributions to PH using specified approximation method
%
% This function scans all service and arrival processes in the network
% structure and converts non-Markovian distributions to Markovian Arrival
% Processes (MAPs) using the specified approximation method.
%
% Input:
%   sn: Network structure from getStruct()
%   options: Solver options structure with fields:
%       - config.nonmkv: Method for conversion ('none', 'bernstein')
%       - config.nonmkvorder: Number of phases for approximation (default 20)
%
% Output:
%   sn: Updated network structure with converted processes
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Get non-Markovian conversion method from options (default 'bernstein')
if isfield(options, 'config') && isfield(options.config, 'nonmkv')
    nonmkvMethod = options.config.nonmkv;
else
    nonmkvMethod = 'bernstein';
end

% If method is 'none', return without any conversion
if strcmpi(nonmkvMethod, 'none')
    return;
end

% Get number of phases from options (default 20)
if isfield(options, 'config') && isfield(options.config, 'nonmkvorder')
    nPhases = options.config.nonmkvorder;
else
    nPhases = 20;
end

% Check if we should preserve deterministic distributions for exact MAP/D/c analysis
if isfield(options, 'config') && isfield(options.config, 'preserveDet')
    preserveDet = options.config.preserveDet;
else
    preserveDet = false;
end

% Markovian ProcessType IDs (no conversion needed)
markovianTypes = [ProcessType.EXP, ProcessType.ERLANG, ProcessType.HYPEREXP, ...
                  ProcessType.PH, ProcessType.APH, ProcessType.MAP, ...
                  ProcessType.COXIAN, ProcessType.COX2, ProcessType.MMPP2, ...
                  ProcessType.IMMEDIATE, ProcessType.DISABLED];

M = sn.nstations;
K = sn.nclasses;

for ist = 1:M
    for r = 1:K
        procType = sn.procid(ist, r);

        % Skip if procType is NaN (e.g., for Transition nodes in SPNs)
        if isnan(procType)
            continue;
        end

        % Skip if already Markovian, disabled, or immediate
        if any(procType == markovianTypes)
            continue;
        end

        % Non-Markovian: need conversion (unless preserveDet for Det)
        distName = ProcessType.toText(procType);
        targetMean = 1 / sn.rates(ist, r);

        % Check if we should skip Det conversion for exact MAP/D/c analysis
        if procType == ProcessType.DET && preserveDet
            % Skip Det - will be handled by exact MAP/D/c solver
            continue;
        end

        % Issue warning for distributions that will be converted
        line_warning(mfilename, ...
            'Distribution %s at station %d class %d is non-Markovian and will be converted to PH (%d phases).\n', ...
            distName, ist, r, nPhases);

        % Get PDF based on distribution type and stored parameters
        origProc = sn.proc{ist}{r};

        switch procType
            case ProcessType.GAMMA
                shape = origProc{1};
                scale = origProc{2};
                pdf_func = @(x) gampdf(x, shape, scale);

            case ProcessType.WEIBULL
                shape_param = origProc{1};  % r
                scale_param = origProc{2};  % alpha
                pdf_func = @(x) wblpdf(x, scale_param, shape_param);

            case ProcessType.LOGNORMAL
                mu = origProc{1};
                sigma = origProc{2};
                pdf_func = @(x) lognpdf(x, mu, sigma);

            case ProcessType.PARETO
                shape_param = origProc{1};  % alpha
                scale_param = origProc{2};  % k (minimum value)
                % Pareto PDF: alpha * k^alpha / x^(alpha+1) for x >= k
                pdf_func = @(x) (x >= scale_param) .* shape_param .* scale_param.^shape_param ./ x.^(shape_param + 1);

            case ProcessType.UNIFORM
                minVal = origProc{1};
                maxVal = origProc{2};
                pdf_func = @(x) (x >= minVal & x <= maxVal) / (maxVal - minVal);

            case ProcessType.DET
                % Deterministic: use Erlang approximation (preserveDet case already handled above)
                MAP = map_erlang(targetMean, nPhases);
                sn = updateSnForMAP(sn, ist, r, MAP, nPhases);
                continue;

            otherwise
                % Generic fallback: Erlang approximation
                MAP = map_erlang(targetMean, nPhases);
                sn = updateSnForMAP(sn, ist, r, MAP, nPhases);
                continue;
        end

        % Apply Bernstein approximation and rescale to target mean
        MAP = map_bernstein(pdf_func, nPhases);
        MAP = map_scale(MAP, targetMean);

        % Update the network structure for the converted MAP
        actualPhases = size(MAP{1}, 1);
        sn = updateSnForMAP(sn, ist, r, MAP, actualPhases);
    end
end

end

function sn = updateSnForMAP(sn, ist, r, MAP, nPhases)
% UPDATESNFORMAP Update all network structure fields for converted MAP
%
% Updates proc, procid, phases, phasessz, phaseshift, mu, phi, pie, nvars, state

% Save old phasessz before updating (needed for state expansion)
oldPhases = sn.phasessz(ist, r);

% Update process representation
sn.proc{ist}{r} = MAP;
% Use PH for Sources (they don't need MAP local variable tracking), MAP for others
if sn.sched(ist) == SchedStrategy.EXT
    sn.procid(ist, r) = ProcessType.PH;
else
    sn.procid(ist, r) = ProcessType.MAP;
end
sn.phases(ist, r) = nPhases;

% Update phasessz and phaseshift (derived from phases)
sn.phasessz(ist, r) = max(nPhases, 1);
% Recompute phaseshift for this station (cumulative sum across classes)
sn.phaseshift(ist, :) = [0, cumsum(sn.phasessz(ist, :))];

% Update mu (rates from -diag(D0))
sn.mu{ist}{r} = -diag(MAP{1});

% Update phi (completion probabilities: sum(D1,2) / -diag(D0))
D0_diag = -diag(MAP{1});
D1_rowsum = sum(MAP{2}, 2);
sn.phi{ist}{r} = D1_rowsum ./ D0_diag;

% Update pie (initial phase distribution)
sn.pie{ist}{r} = map_pie(MAP);

% Update nvars for MAP local variable (skip for Sources - they don't use local vars)
ind = sn.stationToNode(ist);
if sn.sched(ist) ~= SchedStrategy.EXT
    sn.nvars(ind, r) = sn.nvars(ind, r) + 1;
    % Expand state vector (pass old and new phases for proper expansion)
    sn = expandStateForMAP(sn, ind, r, oldPhases, nPhases);
end

% Add PHASE sync event if phases > 1 and not already present
if nPhases > 1
    sn = addPhaseSyncIfNeeded(sn, ind, r);
end
end

function sn = expandStateForMAP(sn, ind, r, oldPhases, newPhases)
% EXPANDSTATEFORMAP Expand state vector for MAP conversion
%
% When converting a non-Markovian distribution to MAP, the state vector
% needs to be expanded in two ways:
% 1. The server portion (space_srv) needs additional columns for extra phases
% 2. The local variable portion (space_var) needs a column for MAP phase memory
%
% State format: [space_buf | space_srv | space_var]
% - space_srv has sum(phasessz) columns total
% - space_var has sum(nvars) columns total

isf = sn.nodeToStateful(ind);
if isf <= 0 || isempty(sn.state) || isempty(sn.state{isf})
    return;
end

ist = sn.nodeToStation(ind);
nRows = size(sn.state{isf}, 1);

% Calculate state vector structure before conversion
% V = number of local variables (already incremented for this class)
V = sum(sn.nvars(ind, :));
V_old = V - 1;  % Before we added the new local var

% K = phases array for this station (already updated for this class)
K = sn.phasessz(ist, :);
sumK_new = sum(K);
K_old = K;
K_old(r) = oldPhases;  % What it was before
sumK_old = sum(K_old);

% Phaseshift tells us where each class's phases start
% For class r, server phases are at positions phaseshift(r)+1 to phaseshift(r)+K(r)
% But phaseshift has already been updated, so compute old positions
Ks_old = [0, cumsum(K_old)];

% Current state dimensions
currentCols = size(sn.state{isf}, 2);

% Calculate buffer size (state columns before space_srv)
% Expected: currentCols = bufSize + sumK_old + V_old
bufSize = currentCols - sumK_old - V_old;
if bufSize < 0
    bufSize = 0;
end

% Extract state portions
if bufSize > 0
    space_buf = sn.state{isf}(:, 1:bufSize);
else
    space_buf = zeros(nRows, 0);
end

if sumK_old > 0
    space_srv = sn.state{isf}(:, bufSize+1:bufSize+sumK_old);
else
    space_srv = zeros(nRows, 0);
end

if V_old > 0
    space_var = sn.state{isf}(:, bufSize+sumK_old+1:end);
else
    space_var = zeros(nRows, 0);
end

% Expand space_srv: insert (newPhases - oldPhases) zeros after class r's position
phasesToAdd = newPhases - oldPhases;
if phasesToAdd > 0
    % Position where class r's phases end (in space_srv)
    insertPos = Ks_old(r) + oldPhases;

    % Insert zeros for the new phases
    space_srv_new = [space_srv(:, 1:insertPos), ...
                     zeros(nRows, phasesToAdd), ...
                     space_srv(:, insertPos+1:end)];
else
    space_srv_new = space_srv;
end

% Expand space_var: add 1 column for the MAP local variable
space_var_new = [space_var, ones(nRows, 1)];

% Reconstruct state
sn.state{isf} = [space_buf, space_srv_new, space_var_new];

% Also update space if it exists
if isfield(sn, 'space') && ~isempty(sn.space) && ~isempty(sn.space{isf})
    nRowsSpace = size(sn.space{isf}, 1);
    currentColsSpace = size(sn.space{isf}, 2);

    % Same calculation for space
    bufSizeSpace = currentColsSpace - sumK_old - V_old;
    if bufSizeSpace < 0
        bufSizeSpace = 0;
    end

    if bufSizeSpace > 0
        space_buf_s = sn.space{isf}(:, 1:bufSizeSpace);
    else
        space_buf_s = zeros(nRowsSpace, 0);
    end

    if sumK_old > 0
        space_srv_s = sn.space{isf}(:, bufSizeSpace+1:bufSizeSpace+sumK_old);
    else
        space_srv_s = zeros(nRowsSpace, 0);
    end

    if V_old > 0
        space_var_s = sn.space{isf}(:, bufSizeSpace+sumK_old+1:end);
    else
        space_var_s = zeros(nRowsSpace, 0);
    end

    if phasesToAdd > 0
        space_srv_s_new = [space_srv_s(:, 1:insertPos), ...
                          zeros(nRowsSpace, phasesToAdd), ...
                          space_srv_s(:, insertPos+1:end)];
    else
        space_srv_s_new = space_srv_s;
    end

    space_var_s_new = [space_var_s, ones(nRowsSpace, 1)];
    sn.space{isf} = [space_buf_s, space_srv_s_new, space_var_s_new];
end
end

function sn = addPhaseSyncIfNeeded(sn, ind, r)
% ADDPHASESYNCIFNEEDED Add a PHASE sync event for converted MAP
%
% When a non-Markovian distribution is converted to a MAP with multiple phases,
% we need to add a PHASE sync event so that phase transitions can occur
% during simulation.

% Check if PHASE sync already exists for this node/class
phaseSyncExists = false;
local = sn.nnodes + 1;

if ~isempty(sn.sync)
    for s = 1:length(sn.sync)
        if ~isempty(sn.sync{s}) && ~isempty(sn.sync{s}.active) && ~isempty(sn.sync{s}.active{1})
            activeEvent = sn.sync{s}.active{1};
            if activeEvent.event == EventType.PHASE && activeEvent.node == ind && activeEvent.class == r
                phaseSyncExists = true;
                break;
            end
        end
    end
end

% Add PHASE sync if not present
if ~phaseSyncExists
    newSync = struct('active', cell(1), 'passive', cell(1));
    newSync.active{1} = Event(EventType.PHASE, ind, r);
    newSync.passive{1} = Event(EventType.LOCAL, local, r, 1.0);
    sn.sync{end+1, 1} = newSync;
end
end
