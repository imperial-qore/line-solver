function [XN,UN,QN,RN,TN,CN,tranSysState,tranSync,sn] = solver_ssa_analyzer_nrm(sn, options)
% SOLVER_SSA_ANALYZER_NRM  Performance indices from SSA/NRM simulation
%   This variant runs the next‑reaction–method SSA on the network SN and
%   returns mean throughput (XN), utilisation (UN), queue length (QN),
%   response time (RN), throughput per station (TN) and class residence
%   time (CN).  Results are based on the empirical state probabilities
%   obtained from SOLVER_SSA_NRM.
%   *** Only stations whose scheduling policy is INF, EXT or PS are
%   supported. Any other policy triggers an error. Cache nodes are not
%   supported in this simplified variant. ***

% Validate scheduling policies ---------------------------------------------------
allowedSched = [SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS];
if any(~arrayfun(@(s) any(s == allowedSched), sn.sched))
    error('solver_ssa_analyzer_nrm:UnsupportedPolicy', ...
        'Only SchedStrategy.INF, SchedStrategy.EXT and SchedStrategy.PS are supported.');
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
if isfield(options, 'config') && isfield(options.config, 'state_space_gen') && ...
        (strcmp(options.config.state_space_gen, 'none') || strcmp(options.config.state_space_gen, 'default'))
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

            case SchedStrategy.PS
                for k = 1:K
                    if ~isempty(PH{ist}{k})
                        UN(ist,k) = pi * depRates(:, (ind-1)*K+k) / rates(ist,k) / S(ist);
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
