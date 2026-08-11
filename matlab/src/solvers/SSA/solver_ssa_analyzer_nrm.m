function [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = solver_ssa_analyzer_nrm(sn, options)
% SOLVER_SSA_ANALYZER_NRM  Performance indices from SSA/NRM simulation
%   This variant runs the next‑reaction–method SSA on the network SN and
%   returns mean throughput (XN), utilisation (UN), queue length (QN),
%   response time (RN), throughput per station (TN) and class residence
%   time (CN).  Results are based on the empirical state probabilities
%   obtained from SOLVER_SSA_NRM.
%   *** Only the scheduling policies listed in ALLOWEDSCHED below are
%   supported. Any other policy triggers an error. Cache nodes are not
%   supported in this simplified variant. ***

% Validate scheduling policies ---------------------------------------------------
% Keep this list in sync with the propensity switch of SOLVER_SSA_NRM and with
% the NRM eligibility gate of SOLVER_SSA_ANALYZER.
allowedSched = [SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS, ...
    SchedStrategy.LPS, SchedStrategy.DPS, SchedStrategy.GPS, ...
    SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO, ...
    SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, ...
    SchedStrategy.HOL, SchedStrategy.SEPT, SchedStrategy.LEPT, ...
    SchedStrategy.LCFSPR, SchedStrategy.PAS, SchedStrategy.POLLING];
if any(~arrayfun(@(s) any(s == allowedSched), sn.sched))
    unsupported = unique(sn.sched(~arrayfun(@(s) any(s == allowedSched), sn.sched)));
    names = strjoin(arrayfun(@(s) SchedStrategy.toText(s), unsupported, 'UniformOutput', false), ', ');
    error('solver_ssa_analyzer_nrm:UnsupportedPolicy', ...
        'The NRM method does not support the scheduling policy: %s.', names);
end

% see _kb/06-solver-catalog.md for rationale (SSA immfeed self-loop)
if isfield(sn,'immfeed') && ~isempty(sn.immfeed) && any(sn.immfeed(:))
    line_warning(mfilename,'SolverSSA(method=nrm) does not model immediate feedback (immfeed); self-loops are treated as class-switching with re-queueing. Use method=''serial'' for immediate feedback.\n');
end

% Shorthands --------------------------------------------------------------------
M    = sn.nstations;               % number of stations
K    = sn.nclasses;                % number of classes
S    = sn.nservers;                % server multiplicities per station
NK   = sn.njobs';                  % population vector (closed chains)
PH   = sn.proc;                    % service‑process MAPs/PHs
sched = sn.sched;                  % scheduling policies

% Pre‑allocate performance vectors/matrices --------------------------------------
XN = NaN(1,K);         % class throughput
UN = NaN(M,K);         % utilisation
QN = NaN(M,K);         % mean jobs at station
RN = NaN(M,K);         % response time
TN = NaN(M,K);         % class throughput at station (departures)
CN = NaN(1,K);         % cycle time per class
% -------------------------------------------------------------------------------
% 1)  Run the SSA/NRM simulator ---------------------------------------------------
% -------------------------------------------------------------------------------
% see _kb/06-solver-catalog.md for rationale (SSA NRM engine dispatch)
useBufferedNrm = true;
if isfield(options, 'config') && isfield(options.config, 'state_space_gen')
    ssg = options.config.state_space_gen;
    if ~(strcmp(ssg, 'none') || strcmp(ssg, 'default'))
        useBufferedNrm = false;
    end
end
if useBufferedNrm
    % Compute only stead-state mean performance indices while running nrm
    [QN, UN, RN, TN, CN, XN, ~, sn] = solver_ssa_nrm(sn, options);
else
    % Use explicit state space generation version

    % -------------------------------------------------------------------------------
    % 1)  Run the SSA/NRM simulator ---------------------------------------------------
    % -------------------------------------------------------------------------------
    [pi,space,depRates,sn] = solver_ssa_nrm_space(sn, options);
    pi=pi(:)';

    % -------------------------------------------------------------------------------
    % 2)  Aggregate performance indices ----------------------------------------------
    % -------------------------------------------------------------------------------
    % Global class throughput ---------------------------------------------------------
    for i=1:M
        for k = 1:K
            refnd = sn.stationToNode(sn.refstat(k));   % reference (sink) stateful index
            XN(k) = pi * depRates(:, (refnd-1)*K+k);
        end
    end

    % Station‑level metrics -----------------------------------------------------------
    rates = sn.rates;
    for ist = 1:M
        ind = sn.stationToNode(ist);   % reference (sink) stateful index
        for k = 1:K
            TN(ist,k) = pi * depRates(:, (ind-1)*K+k);
            QN(ist,k) = pi * space(:, (ind-1)*K + k);
        end

        switch sched(ist)
            case {SchedStrategy.INF, SchedStrategy.EXT}
                % Infinite‑server or external delay: utilisation equals mean jobs
                UN(ist,:) = QN(ist,:);

            case {SchedStrategy.PS, SchedStrategy.FCFS, SchedStrategy.LCFS}
                % see _kb/06-solver-catalog.md for rationale (SSA utilization estimator)
                isCd = ~isempty(sn.cdscaling) && ist <= numel(sn.cdscaling) && ~isempty(sn.cdscaling{ist});
                isJd = ~isempty(sn.jdscaling) && ist <= numel(sn.jdscaling) && ~isempty(sn.jdscaling{ist});
                for k = 1:K
                    if ~isempty(PH{ist}{k})
                        if isCd || isJd
                            % Util = T*S/peak, effective peak = product of the
                            % declared class- and joint-dependence peaks.
                            sdiv = 1;
                            if isCd, sdiv = sdiv * sn.cdscalingpeak(ist,k); end
                            if isJd, sdiv = sdiv * sn.jdscalingpeak(ist,k); end
                        else
                            sdiv = S(ist);
                        end
                        if sdiv > 0
                            UN(ist,k) = pi * depRates(:, (ind-1)*K+k) / rates(ist,k) / sdiv;
                        else
                            UN(ist,k) = 0;
                        end
                    end
                end
        end
    end

    % Response time (Little's law) & cycle time --------------------------------------
    for k = 1:K
        for ist = 1:M
            if TN(ist,k) > 0
                RN(ist,k) = QN(ist,k) / TN(ist,k);
            else
                RN(ist,k) = 0;
            end
        end
        CN(k) = NK(k) / XN(k);
    end
end

% -------------------------------------------------------------------------------
% 3)  Post‑process and clean‑up ----------------------------------------------------
% -------------------------------------------------------------------------------
QN(isnan(QN)) = 0;   UN(isnan(UN)) = 0;   RN(isnan(RN)) = 0;
XN(isnan(XN)) = 0;   TN(isnan(TN)) = 0;   CN(isnan(CN)) = 0;

tranSysState = [];  % transient traces removed in this streamlined version
tranSync     = [];
end
