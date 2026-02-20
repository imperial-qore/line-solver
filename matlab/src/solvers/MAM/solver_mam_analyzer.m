function [QN,UN,RN,TN,CN,XN,runtime,method,totiter,percResults] = solver_mam_analyzer(sn, options)
% [QN,UN,RN,TN,CN,XN,RUNTIME,METHOD] = SOLVER_MAM_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

Tstart = tic;

% Preserve deterministic distributions for exact MAP/D/c analysis
if ~isfield(options, 'config')
    options.config = struct();
end
if ~isfield(options.config, 'preserveDet')
    options.config.preserveDet = true;  % Enable exact MAP/D/c solver
end

% Convert non-Markovian distributions to PH (Det preserved if preserveDet=true)
sn = sn_nonmarkov_toph(sn, options);

% Check if the model is mixed (has both open and closed classes)
isOpen = sn_is_open_model(sn);
isClosed = sn_is_closed_model(sn);
isMixed = isOpen && isClosed;

line_debug('MAM analyzer starting: method=%s, isOpen=%d, isClosed=%d', options.method, isOpen, isClosed);

% Mixed models are supported by the dec.source method

if nargin<2 || isempty(options.config) || ~isfield(options.config,'merge')
    % Preserve existing config fields (like preserveDet) while setting defaults
    if ~isfield(options, 'config') || isempty(options.config)
        options.config = struct();
    end
    if ~isfield(options.config, 'merge')
        options.config.merge = 'super';
    end
    if ~isfield(options.config, 'compress')
        options.config.compress = 'mixture.order1';
    end
    if ~isfield(options.config, 'space_max')
        options.config.space_max = 128;
    end
end

method = options.method;
percResults = []; % Initialize as empty

switch options.method
    case 'dec.mmap'
        % service distribution per class scaled by utilization used as
        % departure process
        line_debug('Using dec.mmap method, calling solver_mam');
        [QN,UN,RN,TN,CN,XN,totiter] = solver_mam(sn, options);
    case {'default', 'dec.source'}
        % Check if network is a valid Fork-Join topology for FJ_codes
        [isFJ, fjInfo] = fj_isfj(sn);

        if isFJ
            % Use FJ_codes for Fork-Join analysis
            if strcmpi(options.method, 'default')
                line_debug('Default method: using FJ_codes for Fork-Join topology\n');
            end
            line_debug('Detected Fork-Join topology, using FJ_codes method');
            [QN,UN,RN,TN,CN,XN,totiter,percResults] = solver_mam_fj(sn, options);
            method = 'qiu';
        else
            % Check if network is a valid BMAP/PH/N/N bufferless retrial queue
            [isRetrial, ~] = qsys_is_retrial(sn);

            % Check if network has reneging (queue abandonment) for MAPMSG
            isReneging = hasRenegingPatience(sn);

            if isRetrial
                % Use BMAP/PH/N/N retrial solver
                if strcmpi(options.method, 'default')
                    line_debug('Default method: using retrial for BMAP/PH/N/N bufferless topology\n');
                end
                line_debug('Detected BMAP/PH/N/N retrial topology, using retrial method');
                [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_retrial(sn, options);
                method = 'retrial';
            elseif isReneging
                % Use MAP/M/s+G reneging solver (MAPMSG)
                if strcmpi(options.method, 'default')
                    line_debug('Default method: using MAPMSG for MAP/M/s+G with reneging\n');
                end
                line_debug('Detected reneging topology, using MAPMSG method');
                [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_retrial(sn, options);
                method = 'reneging';
            else
                % arrival process per chain rescaled by visits at each node
                if strcmpi(options.method, 'default')
                    line_debug('Default method: using dec.source\n');
                end
                line_debug('Using default/dec.source method, calling solver_mam_basic');
                [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic(sn, options);
                method = 'dec.source';
            end
        end
    case 'dec.poisson'
        % analyze the network with Poisson streams
        line_debug('Using dec.poisson method with space_max=1, calling solver_mam_basic');
        options.config.space_max = 1;
        [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_basic(sn, options);
    case 'mna'
        if sn_is_open_model(sn)
            line_debug('Using MNA method for open model, calling solver_mna_open');
            [QN,UN,RN,TN,CN,XN,~,totiter] = solver_mna_open(sn, options);
        elseif sn_is_closed_model(sn)
            line_debug('Using MNA method for closed model, calling solver_mna_closed');
            [QN,UN,RN,TN,CN,XN,~,totiter] = solver_mna_closed(sn, options);
        else
            line_error(mfilename,'The mna method in SolverMAM does not support mixed models.');
        end
    case 'ldqbd'
        % Level-Dependent QBD method for single-class closed networks
        if ~sn_is_closed_model(sn)
            line_error(mfilename,'The ldqbd method requires a closed model.');
        end
        if sn.nclasses ~= 1
            line_error(mfilename,'The ldqbd method requires a single-class model.');
        end
        line_debug('Using LDQBD method for closed model, calling solver_mam_ldqbd');
        [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_ldqbd(sn, options);
    case {'inap', 'inapplus', 'exact'}
        % RCAT-based methods (formerly SolverAG)
        line_debug('Using RCAT method: %s', options.method);
        [QN,UN,RN,TN,CN,XN,totiter] = solver_mam_ag(sn, options);
    otherwise
        line_error(mfilename,'Unknown method: %s', options.method);
end

for i=1:sn.nstations
    switch sn.sched(i)
        case SchedStrategy.EXT
            TN(i,:) = sn.rates(i,:);
    end
end

QN(isnan(QN))=0;
CN(isnan(CN))=0;
RN(isnan(RN))=0;
UN(isnan(UN))=0;
XN(isnan(XN))=0;
TN(isnan(TN))=0;

runtime = toc(Tstart);
end

function isReneging = hasRenegingPatience(sn)
% HASRENEGINGPATIENCE Check if model has reneging/patience configured
%
% Returns true if any queue station has reneging (ImpatienceType.RENEGING)
% configured with a patience distribution.

isReneging = false;

% Check if impatienceClass field exists (ImpatienceType: RENEGING, BALKING)
if ~isfield(sn, 'impatienceClass') || isempty(sn.impatienceClass)
    return;
end

% Check if patienceProc field exists
if ~isfield(sn, 'patienceProc') || isempty(sn.patienceProc)
    return;
end

% Check each station for reneging configuration
for ist = 1:sn.nstations
    for r = 1:sn.nclasses
        if sn.impatienceClass(ist, r) == ImpatienceType.RENEGING
            if ~isempty(sn.patienceProc{ist, r})
                isReneging = true;
                return;
            end
        end
    end
end

end
