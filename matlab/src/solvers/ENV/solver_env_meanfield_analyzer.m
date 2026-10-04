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
% METHOD 'meancov' RUNS HERE TOO, and carries a COVARIANCE beside the mean. The
% propagation is the mean-field one -- same sojourn average, same probOrig
% mixture, same fixed point -- so it lives here rather than in an analyzer of its
% own, which would fork the whole exit/mixture/reseed path to change only the
% object being carried. What is added at each step:
%
%   exit     Cov[Q(T)] = E_T[Cov(Q(t)|t)] + Cov_T[E(Q(t)|t)], the law of total
%            variance over the random sojourn T. The first term is the
%            WITHIN-STAGE covariance, taken from the stage solver when it
%            integrates one (getTranAvgVar, i.e. SolverFLD 'kp' or 'dae') and
%            zero otherwise; the second is the TIMING variance, which the sojourn
%            law contributes whatever the stage solver is. A deterministic
%            sojourn has no timing variance and the exit covariance is C(d).
%   mixture  over the origins h of a switch into e, again by the law of total
%            variance: sum_h p_h (C_h + m_h m_h') - m m'. A convex combination of
%            the C_h ALONE is wrong and silently so: it drops the spread of the
%            per-origin means, which is most of the variance when the stages
%            differ, and the answer still looks like a covariance.
%   reseed   options.config.init_qlen / init_qcov on the stage solver, beside the
%            usual initFromMarginal. Both are in the (station,class) index space
%            ir = (r-1)*M+i, so this analyzer never builds a phase-layout object:
%            the lift into phases belongs to the stage solver that owns the
%            layout (solver_fluid_kp).
%
% A stage solver that reports no covariance and reads no seed still runs, and
% then 'meancov' returns the 'meanfield' means with the environment-induced
% second moment reported beside them (SolverENV.getAvgQLenCov).
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
                        if ~isa(self.solvers{e},'SolverFLD') && ~isa(self.ensemble{e}, 'LayeredNetwork')
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
                % 'meancov' needs the WITHIN-STAGE covariance on top of the
                % means. It costs a second integration of the same stage, which
                % is why only this method asks for it; a stage solver that
                % carries a first moment only says so through
                % supportsTransientVariance and contributes the timing variance
                % alone. QCov.C is (M*K)x(M*K)xnumel(QCov.t), indexed
                % ir = (r-1)*M+i, so post_ needs no phase layout.
                if isMeanCov_(self) && ismethod(self.solvers{e},'supportsTransientVariance') ...
                        && self.solvers{e}.supportsTransientVariance()
                    [tvar, ~, ~, QCovt] = self.solvers{e}.getTranAvgVar();
                    if ~isempty(QCovt)
                        results_e.Tran.Avg.QCov = struct('t', tvar(:), 'C', QCovt);
                    end
                end
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
            meancov = isMeanCov_(self);
            Cexit = cell(1,E); % (M*K)x(M*K) exit covariance of each stage
            for e = 1:E
                Cexit{e} = zeros(M*K, M*K);
            end
            for e=1:E
                if meancov
                    % Independent of the DESTINATION, exactly as Qexit is: the
                    % sojourn weight below reads holdTime{e} alone, because
                    % competing exponentials leave the exit time independent of
                    % where the switch goes.
                    Cexit{e} = stageExitCov_(self, it, e, M, K, detSojourn, dvals(e));
                end
                if isempty(self.results{it,e}) || ~isfield(self.results{it,e}, 'Tran') ...
                        || ~isstruct(self.results{it,e}.Tran.Avg) || ~isfield(self.results{it,e}.Tran.Avg, 'Q')
                    % No transient at all: the stage is constant as far as this
                    % coupling can see, so its exit values are its steady state.
                    [Qc, Uc, Tc] = stageSteady_(self, e, M, K);
                    for h = 1:E
                        Qexit{e,h} = Qc;
                        Uexit{e,h} = Uc;
                        Texit{e,h} = Tc;
                    end
                    continue
                end
                if detSojourn
                    % Deterministic sojourn: exit metrics are the transient at
                    % t=d_e. see _kb/06-solver-catalog.md for rationale
                    % Seeded with the steady state, so a metric the stage
                    % solver carries no trajectory for keeps its constant value
                    % instead of being read as zero. see STAGESTEADY_
                    [Qd, Ud, Td] = stageSteady_(self, e, M, K);
                    for i=1:size(self.results{it,e}.Tran.Avg.Q,1)
                        for r=1:size(self.results{it,e}.Tran.Avg.Q,2)
                            Qir = self.results{it,e}.Tran.Avg.Q{i,r};
                            Uir = self.results{it,e}.Tran.Avg.U{i,r};
                            Tir = self.results{it,e}.Tran.Avg.T{i,r};
                            tref = tranTimeBase_(Qir, Uir, Tir);
                            if ~isempty(tref)
                                Qd(i,r) = tranExitOr_(Qir, tref, [], true, dvals(e), Qd(i,r));
                                Ud(i,r) = tranExitOr_(Uir, tref, [], true, dvals(e), Ud(i,r));
                                Td(i,r) = tranExitOr_(Tir, tref, [], true, dvals(e), Td(i,r));
                            end
                        end
                    end
                    for h = 1:E
                        Qexit{e,h} = Qd; Uexit{e,h} = Ud; Texit{e,h} = Td;
                    end
                    continue
                end
                % Seeded with the steady state, not with zeros: see
                % STAGESTEADY_ for why an absent trajectory is a constant.
                [Qc, Uc, Tc] = stageSteady_(self, e, M, K);
                Qd = Qc; Ud = Uc; Td = Tc;
                % ONE PASS FOR EVERY DESTINATION, as the deterministic-sojourn
                % branch above already does: the weight below reads holdTime{e}
                % and the (i,r) trajectory, and neither depends on h, so the
                % exit averages are the same matrix for all E destinations.
                % Computing them inside the h loop evaluated the sojourn CDF E
                % times per (i,r) per sweep for one answer, and refineForCdf_
                % hands map_cdf a 5000-point grid, so a 24-stage environment
                % paid 24x for every cell of every sweep.
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
                            wir = [0, map_cdf(self.envObj.holdTime{e}, tref(2:end)) - map_cdf(self.envObj.holdTime{e}, tref(1:end-1))]';
                            if ~isnan(wir)
                                Qd(i,r) = tranExitOr_(Qir, tref, wir, false, 0, Qc(i,r));
                                Ud(i,r) = tranExitOr_(Uir, tref, wir, false, 0, Uc(i,r));
                                Td(i,r) = tranExitOr_(Tir, tref, wir, false, 0, Tc(i,r));
                            else
                                % The sojourn weight is undefined here, so
                                % the trajectory cannot be averaged; the
                                % stage's own steady state stands in.
                                Qd(i,r) = Qc(i,r);
                                Ud(i,r) = Uc(i,r);
                                Td(i,r) = Tc(i,r);
                            end
                        end
                    end
                end
                for h = 1:E
                    Qexit{e,h} = Qd; Uexit{e,h} = Ud; Texit{e,h} = Td;
                end
            end

            Qentry = cell(1,E); % average entry queue-length
            Centry = cell(1,E); % entry covariance, 'meancov' only
            for e = 1:E
                % Skip stages where analysis failed (no valid Tran results)
                if isempty(self.results{it,e}) || ~isfield(self.results{it,e}, 'Tran') ...
                        || ~isstruct(self.results{it,e}.Tran.Avg) || ~isfield(self.results{it,e}.Tran.Avg, 'Q')
                    continue
                end
                Qentry{e} = zeros(size(Qexit{e}));
                Sentry = zeros(M*K, M*K); % second moment about zero, 'meancov'
                mentry = zeros(M*K, 1);
                for h=1:E
                    % probability of coming from h to e \times resetFun(Qexit from h to e
                    if self.envObj.probOrig(h,e) > 0
                        mh = self.resetFromMarginal{h,e}(Qexit{h,e});
                        Qentry{e} = Qentry{e} + self.envObj.probOrig(h,e) * mh;
                        if meancov
                            % The reset is an arbitrary map on the means, so the
                            % covariance crosses it by the delta method. Then the
                            % MIXTURE over origins is taken on the SECOND MOMENT,
                            % not on the covariances: a convex combination of the
                            % C_h alone drops the spread of the per-origin means,
                            % which is most of the variance when the stages differ.
                            Rhe = resetJacobian_(self.resetFromMarginal{h,e}, Qexit{h,e}, M, K);
                            Che = Rhe * Cexit{h} * Rhe.';
                            mhv = mh(:);
                            mentry = mentry + self.envObj.probOrig(h,e) * mhv;
                            Sentry = Sentry + self.envObj.probOrig(h,e) * (Che + mhv*mhv.');
                        end
                    end
                end
                if meancov
                    Centry{e} = Sentry - mentry*mentry.';
                    Centry{e} = (Centry{e} + Centry{e}.') / 2; % drop the rounding asymmetry
                end
                if ~isa(self.solvers{e},'SolverFLD') && ~isa(self.ensemble{e}, 'LayeredNetwork')
                    Qentry{e} = self.roundMarginalForDiscreteSolver(Qentry{e}, self.sn{e});
                end
                self.solvers{e}.reset();
                self.ensemble{e}.initFromMarginal(Qentry{e});
                if meancov
                    % The second half of the handoff. initFromMarginal carries the
                    % mean through the MODEL; a covariance has no place on a model
                    % object, so it is handed to the stage SOLVER, in the same
                    % (station,class) index space, together with the mean it
                    % belongs to. init_qlen is not redundant with
                    % initFromMarginal: solver_fluid_kp has its own state layout
                    % and deliberately does not read options.init_sol, so without
                    % it a 'kp' stage would restart empty every iteration.
                    if ~isfield(self.solvers{e}.options, 'config') ...
                            || ~isstruct(self.solvers{e}.options.config)
                        self.solvers{e}.options.config = struct();
                    end
                    self.solvers{e}.options.config.init_qlen = Qentry{e};
                    self.solvers{e}.options.config.init_qcov = Centry{e};
                end
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
                % Seeded with the stage's steady state rather than with
                % zeros: a metric its transient carries no trajectory for is one
                % it holds constant, and calling that zero is what made an open
                % cache model report Q, U and T identically zero. see STAGESTEADY_
                [QExit{e}, UExit{e}, TExit{e}] = stageSteady_(self, e, M, K);
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
                                    QExit{e}(i,r) = tranExitOr_(Qir, tref, [], true, dvals(e), QExit{e}(i,r));
                                    UExit{e}(i,r) = tranExitOr_(Uir, tref, [], true, dvals(e), UExit{e}(i,r));
                                    TExit{e}(i,r) = tranExitOr_(Tir, tref, [], true, dvals(e), TExit{e}(i,r));
                                    continue
                                end
                                [tref, Qir, Uir, Tir] = refineForCdf_(tref, ...
                                    self.envObj.holdTime{e}, Qir, Uir, Tir);
                                w{e} = [0, map_cdf(self.envObj.holdTime{e}, tref(2:end)) - map_cdf(self.envObj.holdTime{e}, tref(1:end-1))]';
                                QExit{e}(i,r) = tranExitOr_(Qir, tref, w{e}, false, 0, QExit{e}(i,r));
                                UExit{e}(i,r) = tranExitOr_(Uir, tref, w{e}, false, 0, UExit{e}(i,r));
                                TExit{e}(i,r) = tranExitOr_(Tir, tref, w{e}, false, 0, TExit{e}(i,r));
                            end
                            % No time base: the seeded steady-state value stands.
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

            % 'meancov' also REPORTS the second moment, mixed over the stages by
            % the same law of total variance the handoff uses. The term
            % sum_e p_e (m_e - m)(m_e - m)' is what the environment itself
            % contributes: two stages with identical within-stage variance but
            % different means still leave the queue length varying, and averaging
            % the per-stage covariances alone would report none of it.
            if isMeanCov_(self)
                Sval = zeros(M*K, M*K);
                mval = zeros(M*K, 1);
                for e = 1:E
                    p = self.envObj.probEnv(e);
                    if ~(p > 0)
                        continue
                    end
                    me = QExit{e}(:);
                    Ce = stageExitCov_(self, it, e, M, K, detSojourn, dvals(e));
                    Sval = Sval + p * (Ce + me*me.');
                    mval = mval + p * me;
                end
                QCov = Sval - mval*mval.';
                QCov = (QCov + QCov.') / 2;
                self.result.Avg.QCov = QCov;
                self.result.Avg.QVar = reshape(max(diag(QCov), 0), M, K);
            end
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

            sn1 = self.sn{1};
            if ~isfield(sn1, 'nodetype')
                return
            end
            cacheNodes = find(sn1.nodetype == NodeType.Cache);
            if isempty(cacheNodes)
                return
            end
            ncaches = numel(cacheNodes);

            % Only fluid stages expose the RMF cache transient this fixed point
            % integrates. Any other stage solver still reports its own cache
            % ratios, so blend THOSE rather than reporting nothing: a silent
            % return leaves the caller with no hit ratio and no reason why.
            allFluid = true;
            for e = 1:E
                if ~isa(self.solvers{e},'SolverFLD')
                    allFluid = false;
                    break
                end
            end
            if ~allFluid
                substituteCacheBlend_(self, E, ref, cacheNodes, K);
                return
            end

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
            ar = cell(1, E); sjc = cell(1, E); exitOcc = cell(1, E);
            % Converged per-stage request rates, kept because the per-list hit
            % split is read off them after the sweep.
            lamc = cell(1, E);
            maxSweep = max(1, self.options.iter_max);
            tol = self.options.iter_tol;
            prevEntryFlat = [];
            for sweep = 1:maxSweep
                for e = 1:E
                    opts_e = self.solvers{e}.options;
                    opts_e.method = 'rmf';
                    opts_e.timespan = tspan{e};
                    % The sojourn average int x dF(t/Lam) / int dF(t/Lam) is
                    % integrated WITH the drift (RMF per-request time mapped to
                    % real time by the cache's request rate Lam), so it does not
                    % depend on the ODE output grid.
                    % see _kb/06-solver-catalog.md for rationale
                    [~, ~, ~, cnodes, ar{e}, ~, lamc{e}, sjall] = ...
                        solver_fld_cacheqn_tran(self.sn{e}, opts_e, entryOcc{e}, self.envObj.holdTime{e});
                    exitOcc{e} = cell(1, ncaches);
                    sjc{e} = cell(1, ncaches);
                    for cc = 1:ncaches
                        cidx = find(cnodes == cacheNodes(cc), 1);
                        if isempty(cidx) || isempty(sjall{cidx})
                            continue
                        end
                        sjc{e}{cc} = sjall{cidx};
                        if ~isempty(sjall{cidx}.xbar) && sjall{cidx}.wtot > 0
                            exitOcc{e}{cc} = sjall{cidx}.xbar;
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
                hitT = zeros(1, K); missT = zeros(1, K); hitLT = [];
                for e = 1:E
                    if numel(sjc{e}) < cc || isempty(sjc{e}{cc})
                        continue
                    end
                    sj = sjc{e}{cc};
                    if ~(sj.wtot > 0) || ~isfinite(sj.wtot) || isempty(sj.xbar)
                        continue
                    end
                    cidx = cc; % cnodes ordering matches cacheNodes across stages
                    pe = self.envObj.probEnv(e);
                    % Sojourn-averaged per-item, per-list occupancy of this stage,
                    % the same average hbar takes, so the list rows sum to it.
                    xbar = sj.xbar;
                    for k = 1:K
                        a = ar{e}(cidx, k);
                        if a <= 0
                            continue
                        end
                        % hit/miss are affine in the occupancy, so their sojourn
                        % averages are those of XBAR.
                        hbar = sj.hitprob(k);
                        mbar = sj.missprob(k);
                        hitT(k)  = hitT(k)  + pe * a * hbar;
                        missT(k) = missT(k) + pe * a * mbar;
                        hbl = local_hit_by_list(xbar, lamc{e}{cidx}, k);
                        if isempty(hbl)
                            continue
                        end
                        if isempty(hitLT)
                            hitLT = zeros(K, numel(hbl));
                        end
                        hitLT(k, :) = hitLT(k, :) + pe * a * hbl;
                    end
                end
                hitprob = NaN(1, K); missprob = NaN(1, K);
                hitproblist = [];
                if ~isempty(hitLT)
                    hitproblist = NaN(K, size(hitLT, 2));
                end
                for k = 1:K
                    tot = hitT(k) + missT(k);
                    if tot > 0
                        hitprob(k)  = hitT(k) / tot;
                        missprob(k) = missT(k) / tot;
                        if ~isempty(hitLT)
                            hitproblist(k, :) = hitLT(k, :) / tot;
                        end
                    end
                end
                node = ref.getNodeByIndex(cacheNodes(cc));
                node.setResultHitProb(hitprob);
                node.setResultMissProb(missprob);
                if ~isempty(hitproblist) && ~all(isnan(hitproblist(:)))
                    node.setResultHitProbList(hitproblist);
                end
            end
end

function bool = isMeanCov_(self)
% Whether this run carries a COVARIANCE beside the mean across a switch. One
% predicate, read by analyze_, post_ and finish_, so the three cannot drift.
bool = isfield(self.options,'method') ...
    && any(strcmpi(self.options.method, {'meancov','env.meancov'}));
end

function R = resetJacobian_(resetFun, Qexit, M, K)
% R = RESETJACOBIAN_(RESETFUN, QEXIT, M, K)
% The Jacobian of a reset policy at the exit mean, so a covariance can cross the
% switch as R*C*R'.
%
% A reset policy is an arbitrary map @(q) -> q on the (station x class) mean
% queue lengths, and there is no general way to push a second moment through one.
% The delta method is the first-order image, which is the order the whole
% mean-field coupling works to. The two NAMED policies are linear and R is then
% EXACT: the identity for 'keep', zero for 'clear'. Linear indexing of the M-by-K
% matrix is column-major, i.e. ir = (r-1)*M+i, the same index space as the
% covariance.
n = M*K;
R = zeros(n, n);
base = resetFun(Qexit);
base = base(:);
if numel(base) ~= n
    return
end
step = 1e-6 * max(1, max(abs(Qexit(:))));
for j = 1:n
    Qp = Qexit;
    Qp(j) = Qp(j) + step;
    pert = resetFun(Qp);
    R(:,j) = (pert(:) - base) / step;
end
end

function C = stageExitCov_(self, it, e, M, K, detSojourn, dval)
% C = STAGEEXITCOV_(SELF, IT, E, M, K, DETSOJOURN, DVAL)
% The covariance of the station-class queue lengths at the instant stage E is
% left, by the law of total variance over the random sojourn T:
%
%   Cov[Q(T)] = E_T[Cov(Q(t)|t)] + Cov_T[E(Q(t)|t)]
%             = sum_n w_n (C(t_n) + m(t_n)m(t_n)')/sum(w) - m_exit m_exit'
%
% with the SAME weights w the exit mean uses, so the two are consistent by
% construction. C(t) is the within-stage covariance the stage solver integrated
% and is absent -- hence zero -- for a stage solver carrying a first moment only;
% the timing term survives regardless. A deterministic sojourn contributes no
% timing variance, so the answer is C(d) alone.
C = zeros(M*K, M*K);
if it < 1 || it > size(self.results,1)
    return
end
res = self.results{it,e};
if isempty(res) || ~isfield(res,'Tran') || ~isstruct(res.Tran.Avg) ...
        || ~isfield(res.Tran.Avg,'Q')
    return
end
% ONE grid for the whole stage. The exit moment multiplies station-class means by
% one another, so they have to be read at the same instants; the per-metric bases
% the exit MEAN uses are each self-contained and need no such alignment.
tref = [];
for i = 1:size(res.Tran.Avg.Q,1)
    for r = 1:size(res.Tran.Avg.Q,2)
        tb = tranTimeBase_(res.Tran.Avg.Q{i,r}, res.Tran.Avg.U{i,r}, res.Tran.Avg.T{i,r});
        if numel(tb) > numel(tref)
            tref = tb(:);
        end
    end
end
if numel(tref) < 2
    return
end
if ~detSojourn
    tref = refineForCdf_(tref(:).', self.envObj.holdTime{e});
    tref = tref(:);
end
nt = numel(tref);
Mtraj = zeros(M*K, nt);
for i = 1:min(M, size(res.Tran.Avg.Q,1))
    for r = 1:min(K, size(res.Tran.Avg.Q,2))
        m = res.Tran.Avg.Q{i,r};
        if isstruct(m) && isfield(m,'metric') && ~isempty(m.metric) ...
                && isfield(m,'t') && numel(m.t) == numel(m.metric) && numel(m.t) > 1
            tm = m.t(:);
            tq = min(max(tref, tm(1)), tm(end));
            Mtraj((r-1)*M + i, :) = interp1(tm, m.metric(:), tq, 'linear').';
        end
    end
end
Ctraj = [];
if isfield(res.Tran.Avg,'QCov')
    Ctraj = interpCov_(res.Tran.Avg.QCov, tref, M*K);
end
if detSojourn
    d = max(tref(1), min(dval, tref(end)));
    if ~isempty(Ctraj)
        Cd = interpCov_(res.Tran.Avg.QCov, d, M*K);
        if ~isempty(Cd)
            C = Cd(:,:,1);
        end
    end
    return
end
trow = tref(:).';
w = [0, map_cdf(self.envObj.holdTime{e}, trow(2:end)) ...
    - map_cdf(self.envObj.holdTime{e}, trow(1:end-1))].';
if numel(w) ~= nt || any(isnan(w)) || ~(sum(w) > 0)
    return
end
sw = sum(w);
mexit = (Mtraj * w) / sw;
S = ((Mtraj .* w.') * Mtraj.') / sw;
if ~isempty(Ctraj)
    S = S + reshape(reshape(Ctraj, (M*K)*(M*K), nt) * w, M*K, M*K) / sw;
end
C = S - mexit*mexit.';
C = (C + C.') / 2;
end

function Ct = interpCov_(qcov, tq, n)
% CT = INTERPCOV_(QCOV, TQ, N)
% A stage's covariance trajectory read at the instants TQ, clamped to its own
% grid rather than extrapolated: linear extrapolation of a covariance can leave
% the positive semidefinite cone, while a convex combination of two members
% stays inside it. Empty when the stage carries no covariance.
Ct = [];
if isempty(qcov) || ~isstruct(qcov) || ~isfield(qcov,'t') || ~isfield(qcov,'C')
    return
end
tv = qcov.t(:);
Craw = qcov.C;
if isempty(tv) || isempty(Craw) || size(Craw,1) ~= n || size(Craw,2) ~= n ...
        || size(Craw,3) ~= numel(tv)
    return
end
tq = tq(:);
if numel(tv) == 1
    Ct = repmat(Craw(:,:,1), 1, 1, numel(tq));
    return
end
tc = min(max(tq, tv(1)), tv(end));
flat = reshape(Craw, n*n, numel(tv)).';
Ct = reshape(interp1(tv, flat, tc, 'linear').', n, n, numel(tq));
end

function hbl = local_hit_by_list(xbar, lam, k)
% Hit probability of class K resolved by cache list, from the mean-field
% occupancy XBAR (flat DDPP state, item-major over lists 0..h) and the
% isolated per-item request rates LAM. List 0 is outside the cache, so the
% hit rows are lists 1..h and sum to the total hit probability of the class.
hbl = [];
if isempty(xbar) || isempty(lam) || k > size(lam, 1)
    return
end
nitems = size(lam, 2);
if nitems <= 0 || mod(numel(xbar), nitems) ~= 0
    return
end
h = numel(xbar) / nitems - 1;
if h < 1
    return
end
wpop = lam(k, :, 1);
wpop(~isfinite(wpop)) = 0;
tot = sum(wpop);
if tot <= 0
    return
end
wpop = wpop(:) / tot;
hbl = zeros(1, h);
for l = 1:h
    hbl(l) = max(0, min(1, wpop' * xbar((1:nitems)' + l * nitems)));
end
end

function substituteCacheBlend_(self, E, ref, cacheNodes, K)
% Cache blend for stages the mean-field cache fixed point cannot integrate.
% Each stage solver has already written its own hit/miss/delayed ratios onto
% its own cache node, so blend those by the per-stage arrival into the cache,
% the same rate weighting the 'avg'/'dec' limits use. The occupancy handoff
% across environment switches is NOT modelled here: every stage is read at its
% own steady state, which is the substitution the warning names.
line_warning(mfilename, ['The mean-field cache fixed point requires fluid stage solvers; ' ...
    'blending each stage''s own cache ratios instead, without the occupancy handoff ' ...
    'across environment switches.']);
nnodes = max(cacheNodes);
acc = SolverENV.newCacheAccumulator(nnodes);
for e = 1:E
    ANn = stageSteadyNodeArv_(self.solvers{e});
    for c = cacheNodes(:)'
        cacheNode = self.ensemble{e}.getNodeByIndex(c);
        acc = SolverENV.accumCacheMetric(acc, c, cacheNode, ...
            SolverENV.cacheArrivalRow(ANn, c, K) * self.envObj.probEnv(e));
    end
end
[cacheHit, cacheMiss, cacheDHit, cacheHitL] = SolverENV.finalizeCacheAccumulator(acc);
for c = cacheNodes(:)'
    node = ref.getNodeByIndex(c);
    if ~isempty(cacheHit{c});  node.setResultHitProb(cacheHit{c});   end
    if ~isempty(cacheMiss{c}); node.setResultMissProb(cacheMiss{c}); end
    if ~isempty(cacheDHit{c}); node.setResultDelayedHitProb(cacheDHit{c}); end
    if ~isempty(cacheHitL{c}); node.setResultHitProbList(cacheHitL{c}); end
end
end

function ANn = stageSteadyNodeArv_(solver)
% The per-node arrival rates of a stage solved in STEADY STATE, which also
% leaves that stage's steady cache ratios on its Cache nodes.
%
% The blend is a steady-state blend, so it must not inherit the stage solver's
% transient horizon: SolverENV leaves a finite options.timespan on a stage
% solver whose user gave one, getAvg refuses that by name, and the blend died
% with "getAvg method does not support the timespan option". The horizon is put
% aside for this solve, as STAGESTEADY_ does, and the transient result is
% dropped first so the node table is the steady one and not a cached transient.
% Both the horizon and the solver's result are put back afterwards.
hasTs = isfield(solver.options, 'timespan');
if hasTs
    tsSaved = solver.options.timespan;
    solver.options.timespan = [Inf, Inf]; % steady: CTMC reads [0,Inf] as a transient
end
resSaved = solver.result;
solver.reset();
try
    [~,~,~,~,ANn] = solver.getAvgNode();
catch me
    if hasTs, solver.options.timespan = tsSaved; end
    solver.result = resSaved;
    rethrow(me);
end
if hasTs, solver.options.timespan = tsSaved; end
solver.result = resSaved;
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

function [v, ok] = tranExit_(m, t, w, detSojourn, d)
% One exit metric: the transient averaged over the holding time, or read at
% the deterministic sojourn d. OK is false when THIS metric has no trajectory,
% which is not the same as its value being zero -- the time base may have come
% from a DIFFERENT metric of the same station, and a present metric is never
% dropped because a sibling is absent. Callers keep their own seed when OK is
% false; see TRANEXITOR_ and STAGESTEADY_.
v = 0; ok = false;
if ~(isstruct(m) && isfield(m,'metric') && ~isempty(m.metric))
    return
end
if numel(m.metric) ~= numel(t)
    return
end
ok = true;
if detSojourn
    v = detEval_(t, m.metric, d);
else
    v = m.metric(:)'*w(:)/sum(w);
end
end

function v = tranExitOr_(m, t, w, detSojourn, d, fallback)
% The transient exit value of M, or FALLBACK when M carries no trajectory. The
% fallback is the stage's own steady state, which for a metric the stage solver
% holds constant IS its exit value.
[v, ok] = tranExit_(m, t, w, detSojourn, d);
if ~ok
    v = fallback;
end
end

function [Qs, Us, Ts] = stageSteady_(self, e, M, K)
% The stage's own STEADY-STATE tables, the exit value of any metric for which
% its transient carries no trajectory.
%
% A METRIC WITH NO TRAJECTORY IS NOT A METRIC WORTH ZERO. A stage solver reports
% a transient only for the quantities its integration actually advances; what it
% leaves out it holds CONSTANT over the stage, and the value at the exit instant
% is then that constant, whatever the sojourn was. Until 2026-09-13 the coupling
% called every such cell 0 and blended the zeros, which on an OPEN CACHE MODEL
% zeroed the whole answer: a Cache is a StatefulNode and a Sink is no station, so
% the only station is the Source, whose throughput neither the fluid nor the CTMC
% stage integrates -- SolverENV then reported Q, U and T identically zero while
% its hit ratios were correct.
%
% READ THE RESULT FIRST, AND ASK WITHOUT A HORIZON. The stage solve that
% produced the transient may already have filled result.Avg; when it has not,
% getAvg has to be asked, and it CANNOT be asked as the stage solver stands --
% SolverENV leaves the stage horizon on options.timespan and getAvg refuses a
% finite one by name. Calling it blind raises, the catch keeps the zeros this
% function exists to replace, and the caller reports an all-zero table. The
% horizon and the untouched transient result are both put back afterwards,
% since the coupling is mid-iteration and reads them.
Qs = zeros(M, K); Us = zeros(M, K); Ts = zeros(M, K);
Qa = []; Ua = []; Ta = [];
solver = self.solvers{e};
% A LAYERED STAGE ANSWERS getAvg IN A DIFFERENT INDEX SPACE than the one this
% fills. (M,K) here is the coupling's BLOCK-DIAGONAL station x class aggregate,
% what getTranHandles lays out and what the exit blend is written in, whereas
% SolverLN.getAvg is getEnsembleAvg: one COLUMN over the LQN's OWN nodes
% (hosts, tasks, entries, activities), carrying no station or class meaning.
% LOCAL_FINITE_BLOCK cannot tell the two apart, so it read LQN node i as
% station i of the first class and pasted an activity's queue length onto an
% off-block cell. Off-block is exactly where the transient overwrites nothing,
% so the stray column survived the probEnv blend and broke the closed-population
% conservation the stages satisfy. getBlockAvg is the steady table in the RIGHT
% space; the shape below is then already (M,K) and the copy is an identity.
if isa(self.ensemble{e}, 'LayeredNetwork')
    [Qb, Ub, Tb] = solver.getBlockAvg();
    Qs = local_finite_block(Qb, M, K);
    Us = local_finite_block(Ub, M, K);
    Ts = local_finite_block(Tb, M, K);
    return
end
res = solver.result;
if ~isempty(res) && isfield(res, 'Avg') && isstruct(res.Avg)
    if isfield(res.Avg, 'Q'), Qa = res.Avg.Q; end
    if isfield(res.Avg, 'U'), Ua = res.Avg.U; end
    if isfield(res.Avg, 'T'), Ta = res.Avg.T; end
end
if isempty(Ta)
    tsSaved = [];
    hasTs = isfield(solver.options, 'timespan');
    if hasTs
        tsSaved = solver.options.timespan;
        solver.options.timespan = [Inf, Inf]; % the default steady request: CTMC reads [0,Inf] as a transient
    end
    try
        [Qa, Ua, ~, Ta] = solver.getAvg();
    catch
        Qa = []; Ua = []; Ta = [];
    end
    if hasTs, solver.options.timespan = tsSaved; end
    solver.result = res;
    if isempty(Ta)
        return
    end
end
Qs = local_finite_block(Qa, M, K);
Us = local_finite_block(Ua, M, K);
Ts = local_finite_block(Ta, M, K);
end

function B = local_finite_block(A, M, K)
% The (M,K) leading block of A with every non-finite entry left at zero: a NaN
% is ABSENT, not a number, and must not propagate into the probEnv-weighted sum.
B = zeros(M, K);
if isempty(A)
    return
end
m = min(M, size(A,1)); k = min(K, size(A,2));
blk = A(1:m, 1:k);
blk(~isfinite(blk)) = 0;
B(1:m, 1:k) = blk;
end
