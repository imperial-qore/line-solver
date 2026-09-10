function varargout = solver_env_meanfield_analyzer(self, phase, it, e)
% SOLVER_ENV_MEANFIELD_ANALYZER  Mean-field (marginal mean queue-length) coupling for SolverENV.
%
% This analyzer implements the default (and historical) SolverENV coupling:
% every stage is solved transiently and reduced to its marginal mean queue
% lengths, which are carried across environment switches via initFromMarginal.
% This mean-field collapse discards the joint state distribution at the handoff;
% see solver_env_statevec_analyzer for the full state-vector alternative. It is
% the default analyzer (options.method other than 'statevec').
%
% Phase dispatch (called by the SolverENV EnsembleSolver hooks):
%   'pre'       pre_(self,it)            -> []
%   'analyze'   analyze_(self,it,e)      -> [results_e, runtime]
%   'post'      post_(self,it)           -> []
%   'finish'    finish_(self)            -> []
%   'converged' converged_(self,it)      -> bool
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

switch phase
    case 'pre'
        pre_(self, it);
    case 'analyze'
        [varargout{1}, varargout{2}] = analyze_(self, it, e);
    case 'post'
        post_(self, it);
    case 'finish'
        finish_(self);
    case 'converged'
        varargout{1} = converged_(self, it);
    otherwise
        line_error(mfilename, sprintf('Unknown meanfield-analyzer phase: %s', phase));
end
end

function bool = converged_(self, it) % convergence test at iteration it
            % BOOL = CONVERGED(IT) % CONVERGENCE TEST AT ITERATION IT
            % Computes max relative absolute difference of queue lengths between iterations
            % Aligned with JAR SolverEnv.converged() implementation

            bool = false;
            if it <= 1
                return
            end

            E = self.getNumberOfModels;
            M = self.ensemble{1}.getNumberOfStations;
            K = self.ensemble{1}.getNumberOfClasses;

            % Check convergence per class (aligned with JAR structure)
            for k = 1:K
                % Build QEntry and QExit matrices (M x E) for this class
                QEntry = zeros(M, E);
                QExit = zeros(M, E);
                for e = 1:E
                    % Skip stages where analysis failed
                    res_curr = self.results{it,e};
                    res_prev = self.results{it-1,e};
                    for i = 1:M
                        if ~isempty(res_curr) && isfield(res_curr, 'Tran') ...
                                && isstruct(res_curr.Tran.Avg) && isfield(res_curr.Tran.Avg, 'Q')
                            Qik_curr = res_curr.Tran.Avg.Q{i,k};
                            if isstruct(Qik_curr) && isfield(Qik_curr, 'metric') && ~isempty(Qik_curr.metric)
                                QExit(i,e) = Qik_curr.metric(1);
                            end
                        end
                        if ~isempty(res_prev) && isfield(res_prev, 'Tran') ...
                                && isstruct(res_prev.Tran.Avg) && isfield(res_prev.Tran.Avg, 'Q')
                            Qik_prev = res_prev.Tran.Avg.Q{i,k};
                            if isstruct(Qik_prev) && isfield(Qik_prev, 'metric') && ~isempty(Qik_prev.metric)
                                QEntry(i,e) = Qik_prev.metric(1);
                            end
                        end
                    end
                end

                % Compute max relative absolute difference using maxpe
                % maxpe computes max(abs(1 - approx./exact)) = max(abs((approx-exact)./exact))
                % This matches JAR's Matrix.maxAbsDiff() implementation
                maxDiff = maxpe(QExit(:), QEntry(:));
                if isempty(maxDiff)
                    maxDiff = 0;  % Handle case where all QEntry values are zero
                end
                if isnan(maxDiff) || isinf(maxDiff)
                    return  % Non-convergence on invalid values
                end
                if maxDiff >= self.options.iter_tol
                    return  % Not converged
                end
            end
            bool = true;
end

function pre_(self, it)
            % PRE(IT)

            if it==1
                for e=list(self)
                    try
                        if isinf(self.getSolver(e).options.timespan(2))
                            [QN,~,~,~] = self.getSolver(e).getAvg();
                        else
                            [QNt,~,~] = self.getSolver(e).getTranAvg();
                            % Handle case where getTranAvg returns NaN for disabled/empty results
                            QN = zeros(size(QNt));
                            for i = 1:size(QNt,1)
                                for k = 1:size(QNt,2)
                                    if isstruct(QNt{i,k}) && isfield(QNt{i,k}, 'metric')
                                        QN(i,k) = QNt{i,k}.metric(end);
                                    else
                                        QN(i,k) = 0; % Default to 0 for NaN/invalid entries
                                    end
                                end
                            end
                        end
                        if ~isa(self.solvers{e}, 'SolverFluid') && ~isa(self.ensemble{e}, 'LayeredNetwork')
                            QN = self.roundMarginalForDiscreteSolver(QN, self.sn{e});
                        end
                        self.ensemble{e}.initFromMarginal(QN);
                    catch
                        % Skip pre-initialization for stages where getAvg fails
                    end
                end
            end
end

function [results_e, runtime] = analyze_(self, it, e)
            % [RESULTS_E, RUNTIME] = ANALYZE(IT, E)
            results_e = struct();
            results_e.('Tran') = struct();
            results_e.Tran.('Avg') = [];
            T0 = tic;
            runtime = toc(T0);
            %% initialize
            try
                % The resolved horizon is set for THIS stage solve and put back
                % after it: a stage solver's own getters branch on whether a
                % horizon is set, and leaving it there changes what the caller's
                % later getEnsembleAvgTables() returns. cdfGrid_ below needs a
                % finite one, and getTranAvg would otherwise resolve it into a
                % local copy nobody else sees. See SolverENV.stageHorizon_.
                tsStage = self.stageHorizon_(e);
                if numel(tsStage) >= 2
                    tsSaved = self.solvers{e}.options.timespan;
                    self.solvers{e}.options.timespan = tsStage;
                    restoreTs = onCleanup(@() self.setStageTimespan_(e, tsSaved)); %#ok<NASGU>
                end
                [Qt,Ut,Tt] = self.ensemble{e}.getTranHandles;
                self.solvers{e}.reset();
                % ASK FOR THE POINTS THE SOJOURN WEIGHT NEEDS, rather than
                % interpolating the integrator's own grid onto them afterwards:
                % interpolation cannot recover resolution the trajectory never
                % had. Honoured by the fluid iteration (options.tranpoints); a
                % stage solver that ignores it still gets refined in post_.
                self.solvers{e}.options.tranpoints = ...
                    cdfGrid_(self.solvers{e}.options.timespan, self.envObj.holdTime{e});
                [QNt,UNt,TNt] = self.solvers{e}.getTranAvg(Qt,Ut,Tt);
                results_e.Tran.Avg.Q = QNt;
                results_e.Tran.Avg.U = UNt;
                results_e.Tran.Avg.T = TNt;
            catch
                % Skip analysis for stages where transient analysis fails
            end
end

function post_(self, it)
            % POST(IT)

            E = self.getNumberOfModels;
            M = self.ensemble{1}.getNumberOfStations;
            K = self.ensemble{1}.getNumberOfClasses;
            [detSojourn, dvals] = sojournConfig_(self, E);
            for e=1:E
                if isempty(self.results{it,e}) || ~isfield(self.results{it,e}, 'Tran') ...
                        || ~isstruct(self.results{it,e}.Tran.Avg) || ~isfield(self.results{it,e}.Tran.Avg, 'Q')
                    for h = 1:E
                        Qexit{e,h} = zeros(M, K);
                        Uexit{e,h} = zeros(M, K);
                        Texit{e,h} = zeros(M, K);
                    end
                    continue
                end
                if detSojourn
                    % Deterministic sojourn: exit metrics are the transient at
                    % t=d_e. see _kb/06-solver-catalog.md for rationale
                    Qd = zeros(M, K); Ud = zeros(M, K); Td = zeros(M, K);
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            Qir = self.results{it,e}.Tran.Avg.Q{i,r};
                            Uir = self.results{it,e}.Tran.Avg.U{i,r};
                            Tir = self.results{it,e}.Tran.Avg.T{i,r};
                            tref = tranTimeBase_(Qir, Uir, Tir);
                            if ~isempty(tref)
                                Qd(i,r) = tranExit_(Qir, tref, [], true, dvals(e));
                                Ud(i,r) = tranExit_(Uir, tref, [], true, dvals(e));
                                Td(i,r) = tranExit_(Tir, tref, [], true, dvals(e));
                            end
                        end
                    end
                    for h = 1:E
                        Qexit{e,h} = Qd; Uexit{e,h} = Ud; Texit{e,h} = Td;
                    end
                    continue
                end
                for h = 1:E
                    Qexit{e,h} = zeros(size(self.results{it,e}.Tran.Avg.Q));
                    Uexit{e,h} = zeros(size(self.results{it,e}.Tran.Avg.U));
                    Texit{e,h} = zeros(size(self.results{it,e}.Tran.Avg.T));
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            Qir = self.results{it,e}.Tran.Avg.Q{i,r};
                            Uir = self.results{it,e}.Tran.Avg.U{i,r};
                            Tir = self.results{it,e}.Tran.Avg.T{i,r};
                            tref = tranTimeBase_(Qir, Uir, Tir);
                            if ~isempty(tref)
                                % The weight lives on the sojourn scale, so the
                                % grid has to as well; see refineForCdf_.
                                [tref, Qir, Uir, Tir] = refineForCdf_(tref, ...
                                    self.envObj.holdTime{e}, Qir, Uir, Tir);
                                % The handoff averages over the SOJOURN, not over the e->h clock:
                                % competing exponentials leave the exit time independent of the
                                % destination. see _kb/06-solver-catalog.md (ENV meanfield)
                                w{e,h} = [0, map_cdf(self.envObj.holdTime{e}, tref(2:end)) - map_cdf(self.envObj.holdTime{e}, tref(1:end-1))]';
                                if ~isnan(w{e,h})
                                    Qexit{e,h}(i,r) = tranExit_(Qir, tref, w{e,h}, false, 0);
                                    Uexit{e,h}(i,r) = tranExit_(Uir, tref, w{e,h}, false, 0);
                                    Texit{e,h}(i,r) = tranExit_(Tir, tref, w{e,h}, false, 0);
                                else
                                    Qexit{e,h}(i,r) = 0;
                                    Uexit{e,h}(i,r) = 0;
                                    Texit{e,h}(i,r) = 0;
                                end
                            else
                                w{e,h} = 0;
                            end
                        end
                    end
                end
            end

            Qentry = cell(1,E); % average entry queue-length
            for e = 1:E
                % Skip stages where analysis failed (no valid Tran results)
                if isempty(self.results{it,e}) || ~isfield(self.results{it,e}, 'Tran') ...
                        || ~isstruct(self.results{it,e}.Tran.Avg) || ~isfield(self.results{it,e}.Tran.Avg, 'Q')
                    continue
                end
                Qentry{e} = zeros(size(Qexit{e}));
                for h=1:E
                    % probability of coming from h to e \times resetFun(Qexit from h to e
                    if self.envObj.probOrig(h,e) > 0
                        Qentry{e} = Qentry{e} + self.envObj.probOrig(h,e) * self.resetFromMarginal{h,e}(Qexit{h,e});
                    end
                end
                if ~isa(self.solvers{e}, 'SolverFluid') && ~isa(self.ensemble{e}, 'LayeredNetwork')
                    Qentry{e} = self.roundMarginalForDiscreteSolver(Qentry{e}, self.sn{e});
                end
                self.solvers{e}.reset();
                self.ensemble{e}.initFromMarginal(Qentry{e});
            end

            % Update transition rates between stages if state-dependent method
            if isfield(self.options, 'method') && strcmp(self.options.method, 'statedep')
                for e = 1:E
                    for h = 1:E
                        if ~isa(self.envObj.env{e,h}, 'Disabled') && ~isempty(self.resetEnvRates{e,h})
                            self.envObj.env{e,h} = self.resetEnvRates{e,h}(...
                                self.envObj.env{e,h}, Qexit{e,h}, Uexit{e,h}, Texit{e,h});
                        end
                    end
                end
                % Reinitialize environment after rate updates
                self.envObj.init();
            end
end

function finish_(self)
            % FINISH()

            it = size(self.results,1); % use last iteration
            E = self.getNumberOfModels;
            M = self.ensemble{1}.getNumberOfStations;
            K = self.ensemble{1}.getNumberOfClasses;
            [detSojourn, dvals] = sojournConfig_(self, E);
            for e=1:E
                QExit{e}=zeros(M, K);
                UExit{e}=zeros(M, K);
                TExit{e}=zeros(M, K);
                if it>0 && ~isempty(self.results{it,e}) && isfield(self.results{it,e}, 'Tran') ...
                        && isstruct(self.results{it,e}.Tran.Avg) && isfield(self.results{it,e}.Tran.Avg, 'Q')
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            Qir = self.results{it,e}.Tran.Avg.Q{i,r};
                            Uir = self.results{it,e}.Tran.Avg.U{i,r};
                            Tir = self.results{it,e}.Tran.Avg.T{i,r};
                            tref = tranTimeBase_(Qir, Uir, Tir);
                            if ~isempty(tref)
                                if detSojourn
                                    % Deterministic sojourn: evaluate at t=d_e.
                                    QExit{e}(i,r) = tranExit_(Qir, tref, [], true, dvals(e));
                                    UExit{e}(i,r) = tranExit_(Uir, tref, [], true, dvals(e));
                                    TExit{e}(i,r) = tranExit_(Tir, tref, [], true, dvals(e));
                                    continue
                                end
                                [tref, Qir, Uir, Tir] = refineForCdf_(tref, ...
                                    self.envObj.holdTime{e}, Qir, Uir, Tir);
                                w{e} = [0, map_cdf(self.envObj.holdTime{e}, tref(2:end)) - map_cdf(self.envObj.holdTime{e}, tref(1:end-1))]';
                                QExit{e}(i,r) = tranExit_(Qir, tref, w{e}, false, 0);
                                UExit{e}(i,r) = tranExit_(Uir, tref, w{e}, false, 0);
                                TExit{e}(i,r) = tranExit_(Tir, tref, w{e}, false, 0);
                            else
                                QExit{e}(i,r) = 0;
                                UExit{e}(i,r) = 0;
                                TExit{e}(i,r) = 0;
                            end
                        end
                    end
                end
                %                 for h = 1:E
                %                     QE{e,h} = zeros(size(self.results{it,e}.Tran.Avg.Q));
                %                     for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                %                         for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                %                             w{e,h} = [0, map_cdf(self.envObj.proc{e}{h}, self.results{it,e}.Tran.Avg.Q{i,r}(2:end,2)) - map_cdf(self.envObj.proc{e}{h}, self.results{it,e}.Tran.Avg.Q{i,r}(1:end-1,2))]';
                %                             if ~isnan(w{e,h})
                %                                 QE{e,h}(i,r) = self.results{it,e}.Tran.Avg.Q{i,r}(:,1)'*w{e,h}/sum(w{e,h});
                %                             else
                %                                 QE{e,h}(i,r) = 0;
                %                             end
                %                         end
                %                     end
                %                 end
            end

            Qval=0*QExit{e};
            Uval=0*UExit{e};
            Tval=0*TExit{e};
            for e=1:E
                Qval = Qval + self.envObj.probEnv(e) * QExit{e}; % to check
                Uval = Uval + self.envObj.probEnv(e) * UExit{e}; % to check
                Tval = Tval + self.envObj.probEnv(e) * TExit{e}; % to check
            end
            self.result.Avg.Q = Qval;
            %    self.result.Avg.R = R;
            %    self.result.Avg.X = X;
            self.result.Avg.U = Uval;
            self.result.Avg.T = Tval;
            %    self.result.Avg.C = C;
            %self.result.runtime = runtime;
            %if self.options.verbose
            %    line_printf('\n');
            %end

            % Cache-hit aggregation across the environment (probEnv-weighted,
            % sojourn-averaged), written onto the reference model's cache nodes.
            % see _kb/06-solver-catalog.md for rationale
            aggregateCacheMeanfield_(self, it, E);
end

function aggregateCacheMeanfield_(self, it, E) %#ok<INUSD>
            % Cache-hit aggregation for fluid inner solvers. Runs a self-
            % contained mean-field fixed point that carries each cache's mean
            % occupancy across environment switches (the cache analog of the
            % queue-length initFromMarginal handoff): stage e is integrated over
            % its sojourn from an entry occupancy that mixes the exit occupancy
            % of its predecessors by probOrig. At convergence the per-class hit
            % throughput is the probEnv-weighted, sojourn-averaged arrival x
            % hit-prob, and the reported hit ratio is hit/(hit+miss).
            ref = self.ensemble{1};
            K = ref.getNumberOfClasses;

            % Only fluid stages expose the RMF cache transient used here.
            for e = 1:E
                if ~isa(self.solvers{e}, 'SolverFluid')
                    return
                end
            end
            sn1 = self.sn{1};
            if ~isfield(sn1, 'nodetype')
                return
            end
            cacheNodes = find(sn1.nodetype == NodeType.Cache);
            if isempty(cacheNodes)
                return
            end
            ncaches = numel(cacheNodes);

            % Finite integration window per stage (fall back to a few mean
            % holding times when the inner solver left the timespan open).
            tspan = cell(1, E);
            for e = 1:E
                ts = self.solvers{e}.options.timespan;
                if isempty(ts) || any(~isfinite(ts))
                    ts = [0, 20 * map_mean(self.envObj.holdTime{e})];
                end
                tspan{e} = ts;
            end

            entryOcc = repmat({cell(1, ncaches)}, 1, E); % entryOcc{e}{c}
            hp = cell(1, E); mp = cell(1, E); ar = cell(1, E);
            tc = cell(1, E); wmass = cell(1, E); exitOcc = cell(1, E);
            maxSweep = max(1, self.options.iter_max);
            tol = self.options.iter_tol;
            prevEntryFlat = [];
            for sweep = 1:maxSweep
                for e = 1:E
                    opts_e = self.solvers{e}.options;
                    opts_e.method = 'rmf';
                    opts_e.timespan = tspan{e};
                    [tc{e}, hp{e}, mp{e}, cnodes, ar{e}, xocc] = ...
                        solver_fld_cacheqn_tran(self.sn{e}, opts_e, entryOcc{e});
                    tt = tc{e}(:)';
                    % Sojourn-weighted exit occupancy per cache for the handoff
                    % (RMF drift mapped from per-request to real time).
                    % see _kb/06-solver-catalog.md for rationale
                    exitOcc{e} = cell(1, ncaches);
                    wmass{e} = cell(1, ncaches);
                    for cc = 1:ncaches
                        cidx = find(cnodes == cacheNodes(cc), 1);
                        if isempty(cidx)
                            continue
                        end
                        Lam = sum(ar{e}(cidx, :));
                        if ~(Lam > 0)
                            continue
                        end
                        treal = tt / Lam;
                        wm = [0, map_cdf(self.envObj.holdTime{e}, treal(2:end)) ...
                                 - map_cdf(self.envObj.holdTime{e}, treal(1:end-1))];
                        wmass{e}{cc} = wm(:);
                        sw = sum(wm);
                        if ~isempty(xocc{cidx}) && sw > 0
                            exitOcc{e}{cc} = xocc{cidx} * wmass{e}{cc} / sw;
                        end
                    end
                end
                % Update each stage's entry occupancy from its predecessors.
                newEntry = repmat({cell(1, ncaches)}, 1, E);
                for e = 1:E
                    for cc = 1:ncaches
                        acc = [];
                        for h = 1:E
                            po = self.envObj.probOrig(h, e);
                            if po > 0 && ~isempty(exitOcc{h}{cc})
                                if isempty(acc)
                                    acc = po * exitOcc{h}{cc};
                                else
                                    acc = acc + po * exitOcc{h}{cc};
                                end
                            end
                        end
                        newEntry{e}{cc} = acc;
                    end
                end
                entryFlat = [];
                for e = 1:E
                    for cc = 1:ncaches
                        entryFlat = [entryFlat; newEntry{e}{cc}(:)]; %#ok<AGROW>
                    end
                end
                entryOcc = newEntry;
                if ~isempty(prevEntryFlat) && numel(prevEntryFlat) == numel(entryFlat)
                    if norm(entryFlat - prevEntryFlat, inf) < tol
                        prevEntryFlat = entryFlat;
                        break
                    end
                end
                prevEntryFlat = entryFlat;
            end

            % Aggregate the converged hit/miss throughputs and write results.
            for cc = 1:ncaches
                hitT = zeros(1, K); missT = zeros(1, K);
                for e = 1:E
                    if numel(wmass{e}) < cc || isempty(wmass{e}{cc})
                        continue
                    end
                    w = wmass{e}{cc};
                    sw = sum(w);
                    if ~(sw > 0) || any(~isfinite(w))
                        continue
                    end
                    cidx = cc; % cnodes ordering matches cacheNodes across stages
                    pe = self.envObj.probEnv(e);
                    for k = 1:K
                        a = ar{e}(cidx, k);
                        if a <= 0
                            continue
                        end
                        hbar = reshape(hp{e}(cidx,k,:), 1, []) * w / sw;
                        mbar = reshape(mp{e}(cidx,k,:), 1, []) * w / sw;
                        hitT(k)  = hitT(k)  + pe * a * hbar;
                        missT(k) = missT(k) + pe * a * mbar;
                    end
                end
                hitprob = NaN(1, K); missprob = NaN(1, K);
                for k = 1:K
                    tot = hitT(k) + missT(k);
                    if tot > 0
                        hitprob(k)  = hitT(k) / tot;
                        missprob(k) = missT(k) / tot;
                    end
                end
                node = ref.getNodeByIndex(cacheNodes(cc));
                node.setResultHitProb(hitprob);
                node.setResultMissProb(missprob);
            end
end

function [detSojourn, dvals] = sojournConfig_(self, E)
% Sojourn option: 'stochastic' (default) averages each stage's transient over
% the random holding-time CDF; 'deterministic' (off by default) evaluates it at
% the fixed mean duration d_e = map_mean(holdTime{e}).
detSojourn = isfield(self.options,'sojourn') && strcmpi(self.options.sojourn,'deterministic');
dvals = zeros(1, E);
if detSojourn
    for e = 1:E
        dvals(e) = map_mean(self.envObj.holdTime{e});
    end
end
end

function v = detEval_(t, metric, d)
% Transient metric at the deterministic sojourn time d (clamped to the grid).
t = t(:); metric = metric(:);
d = max(t(1), min(d, t(end)));
v = interp1(t, metric, d, 'linear');
end

function t = tranTimeBase_(varargin)
% Grid of the first transient metric that carries one, in the order given.
% A station may report one metric and not another -- a Source has no queue
% but does have a throughput -- so the grid cannot be taken from Q alone.
t = [];
for k = 1:numel(varargin)
    m = varargin{k};
    if isstruct(m) && isfield(m,'t') && isfield(m,'metric') && ~isempty(m.t)
        t = m.t;
        return
    end
end
end

function t = cdfGrid_(timespan, holdMap)
% T = CDFGRID_(TIMESPAN, HOLDMAP)
% The instants an exit average should be summed over, given the horizon and the
% sojourn distribution. Empty when the horizon is not finite, i.e. when there is
% no grid to ask for.
%
% The exit metric is the Stieltjes sum sum_k m(t_k)*[F(t_k)-F(t_{k-1})]. Summed
% on the ODE solver's OWN output grid it is an artifact of step placement,
% because that grid is chosen for the HORIZON and not for the sojourn: on
% renv_node_breakdown's DOWN stage 34 of 532 MATLAB points lay below 5*E[S] and
% the exit read 0.876192 against 0.8394 refined.
%
% NINTERP IS 5000 AND THE GRID IS BUILT UNCONDITIONALLY, which is a deliberate
% departure from native python's former rule (500 points, skipped whenever 50 of
% the solver's own already fell inside the support). That escape does not
% converge: with every stage refined the sum runs 0.462260, 0.460580, 0.459704,
% 0.459272, 0.459138, 0.459122 as NINTERP goes 500, 1e3, 2e3, 5e3, 2e4, 5e4, and
% the C++ engine on a 1e5-point uniform grid of its own answers 0.459171. The
% error is first order because the sum reads the RIGHT endpoint and the
% integrand is not small in the tail: an unstable stage (DOWN here serves 0.5
% against arrivals of 0.8) has a queue growing linearly in t, so the mass beyond
% 5*E[S] multiplies a large metric. 5000 points sit 3.3e-4 from the 5e4-point
% value at a twentieth of the cost. All four codebases build the same grid, so
% each then differs only by its own trajectory.
NINTERP = 5000;
t = [];
if isempty(holdMap) || numel(timespan) < 2 ...
        || ~isfinite(timespan(2)) || ~(timespan(2) > timespan(1))
    return
end
t0 = timespan(1);
tend = timespan(2);
meanSojourn = map_mean(holdMap);
if ~(meanSojourn > 0) || ~isfinite(meanSojourn)
    meanSojourn = (tend - t0) / 10;
end
tcdf = min(tend, 5 * meanSojourn);
if tcdf <= t0
    tcdf = tend;
end
ndense = floor(0.9 * NINTERP);
ntail = NINTERP - ndense;
tf = linspace(t0, tcdf, ndense);
if tcdf < tend && ntail > 1
    ttail = linspace(tcdf, tend, ntail + 1);
    tf = [tf, ttail(2:end)];
end
t = tf(:);
end

function [t, varargout] = refineForCdf_(t, holdMap, varargin)
% [T, M1, M2, ...] = REFINEFORCDF_(T, HOLDMAP, M1, M2, ...)
% The FALLBACK path onto the grid cdfGrid_ describes, for a trajectory that was
% not integrated on it.
%
% `analyze_` asks the stage solver for those instants directly
% (`options.tranpoints`), and the fluid iteration honours the request, so this
% returns immediately for a fluid stage: interpolation cannot recover resolution
% a trajectory never had, and asking the integrator costs nothing. A stage
% solver that ignores the request -- a CTMC stage, say -- still arrives here and
% is resampled, which is better than summing on a grid the weight cannot see.
varargout = varargin;
if isempty(t) || numel(t) < 2 || isempty(holdMap)
    return
end
wasrow = isrow(t);
tv = t(:);
tf = cdfGrid_([tv(1), tv(end)], holdMap);
if isempty(tf)
    return
end
% Already integrated on this scale: the dense block of cdfGrid_ runs up to
% tf(ndense), so a trajectory carrying at least that many points below it is
% the requested grid (or finer) and needs nothing done to it.
ndense = floor(0.9 * numel(tf));
if ndense >= 1 && sum(tv <= tf(ndense)) >= ndense
    return
end
for k = 1:numel(varargin)
    m = varargin{k};
    if isstruct(m) && isfield(m,'metric') && ~isempty(m.metric) ...
            && numel(m.metric) == numel(tv)
        mv = interp1(tv, m.metric(:), tf, 'linear');
        m.metric = reshape(mv, size(tf));
        if isfield(m,'t')
            m.t = tf;
        end
        varargout{k} = m;
    end
end
if wasrow
    t = tf.';
else
    t = tf;
end
end

function v = tranExit_(m, t, w, detSojourn, d)
% One exit metric: the transient averaged over the holding time, or read at
% the deterministic sojourn d. Absent metrics are 0; a PRESENT one is never
% dropped because a DIFFERENT metric of the same station is absent.
v = 0;
if ~(isstruct(m) && isfield(m,'metric') && ~isempty(m.metric))
    return
end
if numel(m.metric) ~= numel(t)
    return
end
if detSojourn
    v = detEval_(t, m.metric, d);
else
    v = m.metric(:)'*w(:)/sum(w);
end
end
