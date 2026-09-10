function [Q,U,R,T,C,X,lG,runtime,iter] = solver_ba_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER] = SOLVER_BA_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

T0=tic;
iter = 1;
Q = []; U = [];
R = []; T = [];
C = []; X = [];
lG = NaN;

line_debug('MVA bound analyzer starting: method=%s, nclasses=%d, njobs=%s', options.method, sn.nclasses, mat2str(sn.njobs));

% The STRUCTURAL premises of every bound family -- single-class closed, fully
% closed, single-server, product form -- are asked here (and, for a run that
% enters through SolverBA, in @SolverBA/runAnalyzer with the feature envelope
% beside them), so this run and the report gate SolverBA.supportsModelMethod
% give one answer whichever a caller meets first. Asking at the top is also
% what turned a silent skip into a refusal: the noniterative branches used to
% wrap their whole body in a plain applicability test with no else, so an
% inapplicable model came back as the initial empty matrices and getAvg
% rendered them as a table of zeros over a completed analysis. The premises a
% feature name CAN express (a delay station, an open class, a MAP) stay in
% SolverBA.getMethodFeatureSet instead.
refusal = ba_method_refusal(sn, options.method, options);
if ~isempty(refusal)
    line_error(mfilename, '%s', refusal);
end

switch options.method
    case {'auto.upper','auto.lower'}
        % AUTO composite: evaluate every noniterative bound and keep the
        % tightest side. Feasibility is probed by execution -- a candidate that
        % rejects the model (multiserver, delay station, regime gate) raises and
        % is skipped -- so the list stays correct as families are added.
        line_debug(['Using AUTO composite bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn);
        up = strcmp(options.method,'auto.upper');
        cand = ba_auto_candidates(up);
        Xbest = NaN; best = '';
        for ci = 1:numel(cand)
            oc = options; oc.method = cand{ci};
            try
                [~,~,~,~,~,Xc] = solver_ba_analyzer(sn, oc);
            catch
                continue
            end
            if isempty(Xc) || ~isfinite(Xc(1,1)) || Xc(1,1) <= 0
                continue
            end
            if isnan(Xbest) || (up && Xc(1,1) < Xbest) || (~up && Xc(1,1) > Xbest)
                Xbest = Xc(1,1); best = cand{ci};
            end
        end
        if isnan(Xbest)
            line_error(mfilename,'Method ''%s'' found no feasible bound for this model.', options.method);
        end
        line_debug('AUTO selected %s (X=%g)', best, Xbest);
        X(1,1) = Xbest;
        [Q,U,R,T,C,lG,X] = ba_fill(sn, V, N, Xbest, Zt, D, up);
        runtime=toc(T0);
    case {'lr.upper','lr.lower'}
        % LP-based Linear Reduction bound via mapqn_bnd_lr_pf, distinct from
        % 'qrf.mmi.linear'. see _kb/06-solver-catalog.md for rationale
        line_debug('Using LP linear-reduction bound (%s)', options.method);
        if any(sn.sched == SchedStrategy.INF)
            line_error(mfilename,'Method lr does not support delay (infinite-server) stations.');
        end
        if strcmp(options.method,'lr.lower'), sense = 'min'; else, sense = 'max'; end
        V = sn.visits{1}(:);
        N = sn.njobs(1);
        M = sn.nstations;
        params = struct();
        params.M = M;
        params.N = N;
        params.mu = sn.rates(:,1);
        totV = sum(V);
        params.r = repmat((V(:)')/totV, M, 1);
        params.verbose = false;
        U = zeros(M,1);
        for ti = 1:M
            res = mapqn_bnd_lr_pf(params, ti, sense);
            U(ti,1) = res.objective;
        end
        % Chain throughput implied by the bounded utilizations (U_i = X*V_i/mu_i)
        cand = U(:,1) .* sn.rates(:,1) ./ V;
        cand = cand(isfinite(cand));
        if isempty(cand), X(1,1) = 0; else, X(1,1) = min(cand); end
        T(:,1) = V * X(1,1);
        if strcmp(sense,'max')
            R(:,1) = (1 ./ sn.rates(:,1)) * N;
        else
            R(:,1) = 1 ./ sn.rates(:,1);
        end
        Q(:,1) = T(:,1) .* R(:,1);
        D = V ./ sn.rates(:,1);
        if strcmp(sense,'max'), C(1,1) = N*sum(D); else, C(1,1) = sum(D); end
        runtime=toc(T0);
    case {'mapamva.upper','mapamva.lower'}
        % MAP-AMVA (Casale-Smirni, DSN 2009) via mapqn_bnd_lr_mva: the LP over
        % the exact mean-value balance equations of a closed MAP queueing
        % network. It is the ONLY family here that consumes the correlation
        % between successive services rather than the service mean alone -- its
        % variables are the PER-PHASE queue lengths QN(i,k) and utilizations
        % UN(i,k), so a workload whose burstiness moves the bottleneck between
        % stations is bounded rather than averaged into a renewal process. That
        % is why 'MAP' and 'MMPP2' reach this family's feature set and no other.
        % see _kb/06-solver-catalog.md for rationale
        line_debug('Using MAP-AMVA LP bound (%s)', options.method);
        up = strcmp(options.method,'mapamva.upper');
        if up, sense = 'max'; else, sense = 'min'; end
        [params, perm, Vp, Sp, N] = ba_mapamva_params(sn);
        Mst = params.M;
        % THREE SWEEPS OF THE SAME LP, and the utilization one runs in BOTH
        % senses on purpose. R_i = Q_i/(V_i*X) rises with Q_i and FALLS with X,
        % so an upper bound on the response time pairs Q_i^max with X^min;
        % dividing by X^max on both sides is what would report an upper R below
        % the exact value and break the bracket.
        Umax = zeros(Mst,1); Umin = zeros(Mst,1); Qbnd = zeros(Mst,1);
        for ti = 1:Mst
            rUmax = mapqn_bnd_lr_mva(params, ti, 0, 'max', 'UN');
            rUmin = mapqn_bnd_lr_mva(params, ti, 0, 'min', 'UN');
            rQ    = mapqn_bnd_lr_mva(params, ti, 0, sense, 'QN');
            Umax(ti) = rUmax.objective;
            Umin(ti) = rUmin.objective;
            Qbnd(ti) = rQ.objective;
        end
        % Utilization law U_i = X*V_i*S_i, exact at a single server under ANY
        % service law, so each station turns its own utilization bound into a
        % throughput bound and the tightest of the M survives. A station with no
        % visits or no service time carries no information and is skipped rather
        % than contributing a zero or an Inf.
        load_p = Vp(:) .* Sp(:);
        ok = isfinite(load_p) & load_p > 0;
        if ~any(ok)
            line_error(mfilename, ['Method ''%s'' found no station with both a positive ' ...
                'visit ratio and a positive mean service time.'], options.method);
        end
        Xup = min(Umax(ok) ./ load_p(ok));
        Xlo = max(Umin(ok) ./ load_p(ok));
        if up
            Xb = Xup; Xopp = Xlo;
        else
            Xb = Xlo; Xopp = Xup;
        end
        % Unpermute: the LP orders the stations with the phase-carrying one last.
        inv_perm = zeros(1,Mst); inv_perm(perm) = 1:Mst;
        Vs = Vp(inv_perm); Qs = Qbnd(inv_perm);
        if up
            Us = Umax(inv_perm);
        else
            Us = Umin(inv_perm);
        end
        X(1,1) = Xb;
        T(:,1) = Vs(:) * Xb;
        U(:,1) = Us(:);
        if Xopp > 0
            R(:,1) = Qs(:) ./ (Vs(:) * Xopp);
            % Delay stations are refused above, so the closed-network response
            % time is N/X exactly and the throughput bracket transfers to it
            % directly. Summing the per-station R bounds instead would add M
            % separately-attained maxima and report a looser number.
            C(1,1) = N / Xopp;
        else
            R(:,1) = Inf;
            C(1,1) = Inf;
        end
        Q(:,1) = Qs(:);
        lG = - N*log(Xb);
        runtime=toc(T0);
    case 'aba.upper'
        line_debug('Using ABA upper bound');
        V = sn.visits{1}(:);
        Z = sum(V(sn.sched == SchedStrategy.INF) ./ sn.rates(sn.sched == SchedStrategy.INF));
        D = V(sn.sched ~= SchedStrategy.INF) ./ sn.rates(sn.sched ~= SchedStrategy.INF);
        Dmax = max(D);
        N = sn.nclosedjobs;
        C(1,1) = Z + N * sum(D);
        X(1,1) = min( 1/Dmax, N / (Z + sum(D)));
        T(:,1) = V .* X(1,1);
        R(:,1) = 1 ./ sn.rates * N;
        R(sn.sched == SchedStrategy.INF,1) = 1 ./ sn.rates(sn.sched == SchedStrategy.INF,1);
        Q(:,1) = T(:,1) .* R(:,1);
        U(:,1) = T(:,1) ./ sn.rates;
        U((sn.sched == SchedStrategy.INF),1) = Q((sn.sched == SchedStrategy.INF),1);
        lG = - N*log(X(1,1)); % approx
        runtime=toc(T0);
    case 'aba.lower'
        line_debug('Using ABA lower bound');
        V = sn.visits{1}(:);
        Z = sum(V(sn.sched == SchedStrategy.INF) ./ sn.rates(sn.sched == SchedStrategy.INF));
        D = V(sn.sched ~= SchedStrategy.INF) ./ sn.rates(sn.sched ~= SchedStrategy.INF);
        N = sn.nclosedjobs;
        X(1,1) = N / (Z + N*sum(D));
        C(1,1) = Z + sum(D);
        T(:,1) = V .* X(1,1);
        R(:,1) = 1 ./ sn.rates;
        Q(:,1) = T(:,1) .* R(:,1);
        U(:,1) = T(:,1) ./ sn.rates;
        U((sn.sched == SchedStrategy.INF),1) = Q((sn.sched == SchedStrategy.INF),1);
        lG = - N*log(X(1,1)); % approx
        runtime=toc(T0);
    case 'bjb.upper'
        line_debug('Using BJB upper bound');
        V = sn.visits{1}(:);
        Z = sum(V(sn.sched == SchedStrategy.INF) ./ sn.rates(sn.sched == SchedStrategy.INF));
        D = V(sn.sched ~= SchedStrategy.INF) ./ sn.rates(sn.sched ~= SchedStrategy.INF);
        Dmax = max(D);
        N = sn.nclosedjobs;
        Xaba_upper_1 =  min( 1/Dmax, (N-1) / (Z + sum(D)));
        Xaba_lower_1 =  (N-1) / (Z + (N-1)*sum(D));
        C(1,1) = (Z+sum(D)+max(D)*(N-1-Z*Xaba_lower_1));
        X(1,1) = min(1/Dmax, N / (Z+sum(D)+mean(D)*(N-1-Z*Xaba_upper_1)));
        T(:,1) = V .* X(1,1);
        % RN undefined in the literature so we use ABA upper
        R(:,1) = 1 ./ sn.rates * N;
        %RN = 0*TN;
        %RN(sn.sched ~= SchedStrategy.INF,1) = NaN *  D+ max(D) ./ V(sn.sched ~= SchedStrategy.INF) .* (N-1-Z*Xaba_lower_1) / (sn.nstations - sum(sn.sched == SchedStrategy.INF));
        R(sn.sched == SchedStrategy.INF,1) = 1 ./ sn.rates(sn.sched == SchedStrategy.INF,1);
        Q(:,1) = T(:,1) .* R(:,1);
        U(:,1) = T(:,1) ./ sn.rates;
        U((sn.sched == SchedStrategy.INF),1) = Q((sn.sched == SchedStrategy.INF),1);
        lG = - N*log(X(1,1)); % approx
        runtime=toc(T0);
    case 'bjb.lower'
        line_debug('Using BJB lower bound');
        V = sn.visits{1}(:);
        Z = sum(V(sn.sched == SchedStrategy.INF) ./ sn.rates(sn.sched == SchedStrategy.INF));
        D = V(sn.sched ~= SchedStrategy.INF) ./ sn.rates(sn.sched ~= SchedStrategy.INF);
        Dmax = max(D);
        N = sn.nclosedjobs;
        Xaba_upper_1 =  min( 1/Dmax, (N-1) / (Z + sum(D)));
        Xaba_lower_1 =  (N-1) / (Z + (N-1)*sum(D));
        C(1,1) = (Z+sum(D)+mean(D)*(N-1-Z*Xaba_upper_1));
        X(1,1) = N / (Z+sum(D)+max(D)*(N-1-Z*Xaba_lower_1));
        T(:,1) = V .* X(1,1);
        % RN undefined in the literature so we use ABA lower
        R(:,1) = 1 ./ sn.rates;
        %RN = 0*TN;
        %RN(sn.sched ~= SchedStrategy.INF,1) = NaN * 1 ./ sn.rates(sn.sched ~= SchedStrategy.INF,1) + mean(D) ./ V(sn.sched ~= SchedStrategy.INF) .* (N-1-Z*Xaba_upper_1) / (sn.nstations - sum(sn.sched == SchedStrategy.INF));
        R(sn.sched == SchedStrategy.INF,1) = 1 ./ sn.rates(sn.sched == SchedStrategy.INF,1);
        Q(:,1) = T(:,1) .* R(:,1);
        U(:,1) = T(:,1) ./ sn.rates;
        U((sn.sched == SchedStrategy.INF),1) = Q((sn.sched == SchedStrategy.INF),1);
        lG = - N*log(X(1,1)); % approx
        runtime=toc(T0);
    case 'pb.upper'
        line_debug('Using PB upper bound');
        V = sn.visits{1}(:);
        Z = sum(V(sn.sched == SchedStrategy.INF) ./ sn.rates(sn.sched == SchedStrategy.INF));
        D = V(sn.sched ~= SchedStrategy.INF) ./ sn.rates(sn.sched ~= SchedStrategy.INF);
        Dmax = max(D);
        N = sn.nclosedjobs;
        Xaba_upper_1 =  min( 1/Dmax, (N-1) / (Z + sum(D)));
        Xaba_lower_1 =  (N-1) / (Z + (N-1)*sum(D));
        Dpb2 = sum(D.^2)/sum(D);
        DpbN = sum(D.^N)/sum(D.^(N-1));
        C(1,1) = (Z+sum(D)+DpbN*(N-1-Z*Xaba_lower_1));
        X(1,1) = min(1/Dmax, N / (Z+sum(D)+Dpb2*(N-1-Z*Xaba_upper_1)));
        T(:,1) = V .* X(1,1);
        % RN undefined in the literature so we use ABA upper
        R(:,1) = 1 ./ sn.rates * N;
        %RN = 0*TN;
        %RN(sn.sched ~= SchedStrategy.INF,1) = NaN * 1 ./ sn.rates(sn.sched ~= SchedStrategy.INF,1) + (D.^N/sum(D.^(N-1))) ./ V(sn.sched ~= SchedStrategy.INF)  * (N-1-Z*Xaba_upper_1);
        R(sn.sched == SchedStrategy.INF,1) = 1 ./ sn.rates(sn.sched == SchedStrategy.INF,1);
        Q(:,1) = T(:,1) .* R(:,1);
        U(:,1) = T(:,1) ./ sn.rates;
        U((sn.sched == SchedStrategy.INF),1) = Q((sn.sched == SchedStrategy.INF),1);
        lG = - N*log(X(1,1)); % approx
        runtime=toc(T0);
    case 'pb.lower'
        line_debug('Using PB lower bound');
        V = sn.visits{1}(:);
        Z = sum(V(sn.sched == SchedStrategy.INF) ./ sn.rates(sn.sched == SchedStrategy.INF));
        D = V(sn.sched ~= SchedStrategy.INF) ./ sn.rates(sn.sched ~= SchedStrategy.INF);
        Dmax = max(D);
        N = sn.nclosedjobs;
        Xaba_upper_1 =  min( 1/Dmax, (N-1) / (Z + sum(D)));
        Xaba_lower_1 =  (N-1) / (Z + (N-1)*sum(D));
        Dpb2 = sum(D.^2)/sum(D);
        DpbN = sum(D.^N)/sum(D.^(N-1));
        C(1,1) = (Z+sum(D)+Dpb2*(N-1-Z*Xaba_upper_1));
        X(1,1) = N / (Z+sum(D)+DpbN*(N-1-Z*Xaba_lower_1));
        T(:,1) = V .* X(1,1);
        % RN undefined in the literature so we use ABA lower
        R(:,1) = 1 ./ sn.rates;
        %RN = 0*TN;
        %RN(sn.sched ~= SchedStrategy.INF,1) = NaN *  1 ./ sn.rates(sn.sched ~= SchedStrategy.INF,1) + (D.^2/sum(D)) ./ V(sn.sched ~= SchedStrategy.INF)  * (N-1-Z*Xaba_upper_1);
        R(sn.sched == SchedStrategy.INF,1) = 1 ./ sn.rates(sn.sched == SchedStrategy.INF,1);
        Q(:,1) = T(:,1) .* R(:,1);
        U(:,1) = T(:,1) ./ sn.rates;
        U((sn.sched == SchedStrategy.INF),1) = Q((sn.sched == SchedStrategy.INF),1);
        lG = - N*log(X(1,1)); % approx
        runtime=toc(T0);    
    case 'sb.upper'
        line_debug('Using SB upper bound');
        if any(sn.sched == SchedStrategy.INF)
            line_error(mfilename,'Unsupported method for a model with infinite-server stations.');
        end
        V = sn.visits{1}(:);
        D = V(sn.sched ~= SchedStrategy.INF) ./ sn.rates(sn.sched ~= SchedStrategy.INF);
        Z = sum(V(sn.sched == SchedStrategy.INF) ./ sn.rates(sn.sched == SchedStrategy.INF));
        N = sn.nclosedjobs;
        A3 = sum(D.^3);
        A2 = sum(D.^2);
        A1 = sum(D);
         Dmax = max(D);
        % Harel UB(n) is defined for n <= N only; the level-3 coefficient
        % is not a bound at N < 3, so fall back to UB(2) (Dallery), which
        % is exact at N <= 2 since UB(N) = TH(N).
        if N >= 3
            cub3 = (A1*A2+A3)/(A1^2+A2);
        else
            cub3 = A2/A1;
        end
        C(1,1) = Z+A1+(N-1)*cub3;
        X(1,1) = min([1/Dmax,N / (Z+A1+(N-1)*cub3)]);
        T(:,1) = V .* X(1,1);
        % RN is undefined in the literature for this bound, so it carries
        % the ABA PESSIMISTIC residence: an upper method must publish an
        % upper-consistent R, else Q = T.*R lands below exact (E[n_i] <=
        % N*U_i makes N/mu_i a valid residence bound)
        R(:,1) = 1 ./ sn.rates * N;
        %RN = 0*TN;
        %RN(sn.sched ~= SchedStrategy.INF,1) = NaN *  1 ./ sn.rates(sn.sched ~= SchedStrategy.INF,1) + (D.^2/sum(D)) ./ V(sn.sched ~= SchedStrategy.INF)  * (N-1-Z*Xaba_upper_1);
        R(sn.sched == SchedStrategy.INF,1) = 1 ./ sn.rates(sn.sched == SchedStrategy.INF,1);
        Q(:,1) = T(:,1) .* R(:,1);
        U(:,1) = T(:,1) ./ sn.rates;
        U((sn.sched == SchedStrategy.INF),1) = Q((sn.sched == SchedStrategy.INF),1);
        lG = - N*log(X(1,1)); % approx
        runtime=toc(T0);        
    case 'sb.lower'
        line_debug('Using SB lower bound');
        if any(sn.sched == SchedStrategy.INF)
            line_error(mfilename,'Unsupported method for a model with infinite-server stations.');
        end
        V = sn.visits{1}(:);
        D = V(sn.sched ~= SchedStrategy.INF) ./ sn.rates(sn.sched ~= SchedStrategy.INF);
        Z = sum(V(sn.sched == SchedStrategy.INF) ./ sn.rates(sn.sched == SchedStrategy.INF));
        N = sn.nclosedjobs;
        AN = sum(D.^N);
        A1 = sum(D);
        % (N-1)*(AN/A1)^(1/(N-1)) -> 0 as N -> 1, leaving the exact
        % single-job cycle time; evaluated directly it divides by zero
        if N == 1
            cterm = 0;
        else
            cterm = (N-1)*(AN/A1)^(1/(N-1));
        end
        C(1,1) = Z+A1+cterm;
        X(1,1) = N / (Z+A1+cterm);
        T(:,1) = V .* X(1,1);
        % RN undefined in the literature so we use ABA lower
        R(:,1) = 1 ./ sn.rates;
        %RN = 0*TN;
        %RN(sn.sched ~= SchedStrategy.INF,1) = NaN *  1 ./ sn.rates(sn.sched ~= SchedStrategy.INF,1) + (D.^2/sum(D)) ./ V(sn.sched ~= SchedStrategy.INF)  * (N-1-Z*Xaba_upper_1);
        R(sn.sched == SchedStrategy.INF,1) = 1 ./ sn.rates(sn.sched == SchedStrategy.INF,1);
        Q(:,1) = T(:,1) .* R(:,1);
        U(:,1) = T(:,1) ./ sn.rates;
        U((sn.sched == SchedStrategy.INF),1) = Q((sn.sched == SchedStrategy.INF),1);
        lG = - N*log(X(1,1)); % approx
        runtime=toc(T0);
    case 'gb.upper'
        line_debug('Using GB upper bound');
        V = sn.visits{1}(:);
        Z = sum(V(sn.sched == SchedStrategy.INF) ./ sn.rates(sn.sched == SchedStrategy.INF));
        D = V(sn.sched ~= SchedStrategy.INF) ./ sn.rates(sn.sched ~= SchedStrategy.INF);
        N = sn.nclosedjobs;
        Dmax = max(D);
        X(1,1) = min(1/Dmax, pfqn_xzgsbup(D,N,Z));
        C(1,1) = N / pfqn_xzgsblow(D,N,Z);
        T(:,1) = V .* X(1,1);
        XNlow = pfqn_xzgsblow(D,N,Z);
        k = 0;
        for i=1:sn.nstations
            if sn.sched(i) == SchedStrategy.INF
                R(i,1) = 1 / sn.rates(i);
                % T(i,1), not X: the delay is visited V(i) times per cycle,
                % and dropping that factor lets Q exceed the population
                Q(i,1) = T(i,1) * R(i,1);
            else
                k = k + 1;
                Q(i,1) = pfqn_qzgbup(D,N,Z,k);
                R(i,1) = Q(i,1) / XNlow / V(i) ;
            end
        end
        R(sn.sched == SchedStrategy.INF,1) = 1 ./ sn.rates(sn.sched == SchedStrategy.INF,1);
        U(:,1) = T(:,1) ./ sn.rates;
        U((sn.sched == SchedStrategy.INF),1) = Q((sn.sched == SchedStrategy.INF),1);
        lG = - N*log(X(1,1)); % approx
        runtime=toc(T0);
    case 'gb.lower'
        line_debug('Using GB lower bound');
        V = sn.visits{1}(:);
        Z = sum(V(sn.sched == SchedStrategy.INF) ./ sn.rates(sn.sched == SchedStrategy.INF));
        D = V(sn.sched ~= SchedStrategy.INF) ./ sn.rates(sn.sched ~= SchedStrategy.INF);
        N = sn.nclosedjobs;
        X(1,1) = pfqn_xzgsblow(D,N,Z);
        C(1,1) = N / pfqn_xzgsbup(D,N,Z);
        T(:,1) = V .* X(1,1);
        XNup = pfqn_xzgsbup(D,N,Z);
        k = 0;
        for i=1:sn.nstations
            if sn.sched(i) == SchedStrategy.INF
                R(i,1) = 1 / sn.rates(i);
                % T(i,1), not X: the delay is visited V(i) times per cycle,
                % and dropping that factor lets Q exceed the population
                Q(i,1) = T(i,1) * R(i,1);
            else
                k = k + 1;
                Q(i,1) = pfqn_qzgblow(D,N,Z,k);
                R(i,1) = Q(i,1) / XNup / V(i) ;
            end
        end
        U(:,1) = T(:,1) ./ sn.rates;
        U((sn.sched == SchedStrategy.INF),1) = Q((sn.sched == SchedStrategy.INF),1);
        lG = - N*log(X(1,1)); % approx
        runtime=toc(T0);
    case {'harel.upper','harel.lower'}
        % Sharp bounds of Harel-Namn-Sturm, distinct from the 'sb' family of
        % the same paper: they extrapolate from the EXACT normalizing constant
        % at populations n <= 7 instead of using the first three power sums.
        line_debug(['Using Harel-Namn-Sturm bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn);
        if Zt > 0
            line_error(mfilename,'Method ''%s'' does not support think times (infinite-server stations).', options.method);
        end
        maxUB = min(N, 7);
        [LB,UB,TH] = pfqn_harel_bounds(D, N, 0, maxUB);
        up = strcmp(options.method,'harel.upper');
        if up
            if maxUB >= 2
                Xb = min(1/max(D), UB(maxUB));
            else
                Xb = min(1/max(D), TH(1));
            end
        else
            Xb = LB;
        end
        [Q,U,R,T,C,lG,X] = ba_fill(sn, V, N, Xb, Zt, D, up);
        runtime=toc(T0);
    case {'mwba.upper','mwba.lower'}
        line_debug(['Using Majumdar-Woodside robust box bound: ' options.method]);
        [Lchain,STchain,Vchain,alpha,Nchain] = sn_get_demands_chain(sn);
        isdelay = (sn.sched == SchedStrategy.INF);
        Zc = sum(Lchain(isdelay,:),1);           % think time per chain
        Vq = Vchain(~isdelay,:);                 % queueing stations
        Sq = STchain(~isdelay,:);
        % map station scheduling to MW discipline codes
        % 0=FIFO, 1=PS, 2=non-preemptive priority, 3=preemptive priority
        qstat = find(~isdelay);
        schedq = zeros(numel(qstat),1);
        for kk=1:numel(qstat)
            schedq(kk) = mwrbb_disc_code(sn.sched(qstat(kk)));
        end
        % per-chain priority (lower value = higher priority); use the
        % reference class priority of each chain
        prioc = zeros(1,sn.nchains);
        for c=1:sn.nchains
            inch = sn.inchain{c};
            if isfield(sn,'refclass') && ~isempty(sn.refclass) && sn.refclass(c)>0
                prioc(c) = sn.classprio(sn.refclass(c));
            else
                prioc(c) = min(sn.classprio(inch));
            end
        end
        [Xlo,Xup,Wlo] = pfqn_mwrbb(Vq,Sq,Nchain,Zc,schedq,prioc);
        if strcmp(options.method,'mwba.upper')
            Xchain = Xup;
        else
            Xchain = Xlo;
        end
        Tchain = zeros(sn.nstations,sn.nchains);
        Uchain = zeros(sn.nstations,sn.nchains);
        Qchain = zeros(sn.nstations,sn.nchains);
        for c=1:sn.nchains
            for i=1:sn.nstations
                Tchain(i,c) = Xchain(c) * Vchain(i,c);
                Uchain(i,c) = Xchain(c) * Lchain(i,c);   % utilization law
                % Q is filled below by ba_chain_qfill: neither the
                % no-contention residence Sq nor the Theorem-1 residence
                % Wlo yields a queue length on the declared side.
            end
        end
        Qchain = ba_chain_qfill(Uchain, Xchain, Lchain, Nchain, isdelay, ...
            strcmp(options.method,'mwba.upper'));
        [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, Qchain, Uchain, [], Tchain, [], Xchain);
        lG = NaN;
        runtime=toc(T0);
    case {'pbh.upper','pbh.lower'}
        line_debug(['Using PBH (Eager-Sevcik) bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn);
        [Xlo,Xhi] = pfqn_pbh(D, N, Zt, ba_level(options));
        up = strcmp(options.method,'pbh.upper');
        [Q,U,R,T,C,lG,X] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'cbh.upper','cbh.lower'}
        line_debug(['Using CBH (Dowdy) bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn);
        [Xlo,Xhi] = pfqn_cbh(D, N, Zt, ba_level(options));
        up = strcmp(options.method,'cbh.upper');
        [Q,U,R,T,C,lG,X] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'pbk.upper','pbk.lower'}
        line_debug(['Using PB(k) bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn);
        [Xlo,Xhi] = pfqn_pbk(D, N, Zt, ba_level(options));
        up = strcmp(options.method,'pbk.upper');
        [Q,U,R,T,C,lG,X] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'bjbk.upper','bjbk.lower'}
        line_debug(['Using BJB(k) bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn);
        [Xlo,Xhi] = pfqn_bjbk(D, N, Zt, ba_level(options));
        up = strcmp(options.method,'bjbk.upper');
        [Q,U,R,T,C,lG,X] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'ssd.upper','ssd.lower'}
        line_debug(['Using SSD (Suri-Dallery) multiserver bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn);
        cvec = sn.nservers(sn.sched ~= SchedStrategy.INF);
        [Xlo,Xhi] = pfqn_ssd(D, N, Zt, cvec(:));
        up = strcmp(options.method,'ssd.upper');
        [Q,U,R,T,C,lG,X] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'cub.upper','mbjb.lower'}
        % cub.upper: Kerola composite upper bound (upper-only in the paper).
        % mbjb.lower: the multiclass Balanced Job Bounds lower bound (Kerola
        % eq. 10) that seeds it; the natural multiclass counterpart of bjb.
        line_debug(['Using multiclass composite/BJB bound: ' options.method]);
        [Lchain,STchain,Vchain,alpha,Nchain] = sn_get_demands_chain(sn);
        isdelay = (sn.sched == SchedStrategy.INF);
        Zc = sum(Lchain(isdelay,:),1);
        Lq = Lchain(~isdelay,:);
        [Xub,Xlb] = pfqn_mcub(Lq, Nchain, Zc);
        up = strcmp(options.method,'cub.upper');
        Xchain = up*Xub(:)' + (~up)*Xlb(:)';
        Tchain = zeros(sn.nstations,sn.nchains);
        Uchain = zeros(sn.nstations,sn.nchains);
        Qchain = zeros(sn.nstations,sn.nchains);
        for c=1:sn.nchains
            Tchain(:,c) = Xchain(c) * Vchain(:,c);
            Uchain(:,c) = Xchain(c) * Lchain(:,c);   % utilization law
        end
        Qchain = ba_chain_qfill(Uchain, Xchain, Lchain, Nchain, isdelay, up);
        [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, Qchain, Uchain, [], Tchain, [], Xchain);
        lG = NaN;
        runtime=toc(T0);
    case {'looping.upper','looping.lower'}
        % Eager Looping: the multiclass bracket that initializes the
        % multiple-class PBH. Pessimistic side eq. (2.23), optimistic side
        % from the response-time lower bound of eq. (2.24).
        line_debug(['Using Eager Looping bound: ' options.method]);
        [Lchain,STchain,Vchain,alpha,Nchain] = sn_get_demands_chain(sn);
        isdelay = (sn.sched == SchedStrategy.INF);
        Zc = sum(Lchain(isdelay,:),1);
        Lq = Lchain(~isdelay,:);
        [Xlo,Xup] = pfqn_looping(Lq, Nchain, Zc);
        up = strcmp(options.method,'looping.upper');
        Xchain = up*Xup(:)' + (~up)*Xlo(:)';
        Tchain = zeros(sn.nstations,sn.nchains);
        Uchain = zeros(sn.nstations,sn.nchains);
        for c=1:sn.nchains
            Tchain(:,c) = Xchain(c) * Vchain(:,c);
            Uchain(:,c) = Xchain(c) * Lchain(:,c);   % utilization law
        end
        Qchain = ba_chain_qfill(Uchain, Xchain, Lchain, Nchain, isdelay, up);
        [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, Qchain, Uchain, [], Tchain, [], Xchain);
        lG = NaN;
        runtime=toc(T0);
    case {'sib.upper','sib.lower'}
        line_debug(['Using SIB (Srinivasan) bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn);
        if Zt > 0
            line_error(mfilename,'Method ''%s'' supports Z=0 (no delay station) only; delay needs the SIB Section-3.2 extension.', options.method);
        end
        [Xlo,Xhi] = pfqn_sib(D, N, 0, ba_level(options));
        up = strcmp(options.method,'sib.upper');
        [Q,U,R,T,C,lG,X] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'scb.upper','scb.lower'}
        % Single-class bounds of Dowdy et al. (1992). THE BRACKETED OBJECT IS
        % NOT THIS MODEL: scb brackets the multiclass system that this
        % single-class model aggregates, so scb.lower is the EXACT single-class
        % throughput and scb.upper adds the demand-free Expression-(3) gap.
        % That is why scb is absent from ba_auto_candidates -- mixing it with
        % families that bracket this model's own solution would compare two
        % different quantities.
        line_debug(['Using SCB (Dowdy et al. 1992) single-class bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn);
        if Zt > 0
            line_error(mfilename,'Method ''%s'' supports Z=0 (no delay station) only; Theorem 3 rests on the delay-free balanced-network throughput.', options.method);
        end
        [Xlo,Xhi] = pfqn_scb(D, N);
        up = strcmp(options.method,'scb.upper');
        [Q,U,R,T,C,lG,X] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'ldbcmp.lower'}
        line_debug(['Using LD-BCMP (Anselmi-Cremonesi) lower throughput bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn);
        % Fixed-rate / delay parameterization (Heffes c=0); LD stations are
        % treated as fixed-rate at their limiting demand.
        [Xlo,~,Qhat] = pfqn_ldbcmp(D, N, Zt, zeros(numel(D),1));
        if isnan(Xlo)
            line_error(mfilename,'Method ''%s'' requires the asymptotic regime N >= Qhat (Qhat=%.4f > N=%d).', options.method, Qhat, N);
        end
        [Q,U,R,T,C,lG,X] = ba_fill(sn, V, N, Xlo, Zt, D, false);
        runtime=toc(T0);
end
end

% Single-class demand extraction for the hierarchical bound cases: returns
% visit vector V, aggregate think time Zt, per-queue demand vector D and the
% closed population N. It no longer tests that the model IS single-class
% closed: BA_METHOD_REFUSAL owns that rule for every caller, and a second
% copy here is what would let the report gate and this run disagree.
function [V,Zt,D,N] = ba_sc_demands(sn)
V = sn.visits{1}(:);
isinf_ = (sn.sched == SchedStrategy.INF);
Zt = sum(V(isinf_) ./ sn.rates(isinf_));
D  = V(~isinf_) ./ sn.rates(~isinf_);
N  = sn.nclosedjobs;
end

% Candidate list for the AUTO composite. Noniterative families only: the
% level-parameterized hierarchies (pbh/cbh/pbk/bjbk/sib) and the LP reductions
% are excluded because their cost is not O(K) and their accuracy is a user
% choice, not a fixed property.
% Per-chain queue lengths from a chain-throughput bound, on the declared side.
% Lower: Q_ic >= U_ic, since the station holds a class-c job whenever it serves
% one. Upper: E[n_ic] <= N_c*P(station busy) = N_c*min(1,sum_c U_ic), and a
% delay station queues nothing, so there Q_ic = X_c*L_ic exactly.
function Qchain = ba_chain_qfill(Uchain, Xchain, Lchain, Nchain, isdelay, isUpper)
[M,C] = size(Uchain);
Qchain = zeros(M,C);
Utot = min(1, sum(Uchain,2));
for c=1:C
    for i=1:M
        if isdelay(i)
            Qchain(i,c) = Xchain(c) * Lchain(i,c);
        elseif isUpper
            Qchain(i,c) = Nchain(c) * Utot(i);
        else
            Qchain(i,c) = Uchain(i,c);
        end
    end
end
end

function cand = ba_auto_candidates(isUpper)
if isUpper
    cand = {'aba.upper','bjb.upper','pb.upper','gb.upper','sb.upper', ...
        'mwba.upper','ssd.upper','cub.upper'};
else
    cand = {'aba.lower','bjb.lower','pb.lower','gb.lower','sb.lower', ...
        'mwba.lower','ssd.lower','mbjb.lower','ldbcmp.lower'};
end
end

function lvl = ba_level(options)
lvl = 2;
if isfield(options,'level') && ~isempty(options.level)
    lvl = options.level;
end
end

% Fill per-station [Q,U,R,T,C] from a scalar chain-throughput bound X, using
% the same optimistic (isUpper) / pessimistic residence construction as the
% ABA bound so conventions match across the bound families.
function [Q,U,R,T,C,lG,Xout] = ba_fill(sn, V, N, X, Zt, D, isUpper)
isinf_ = (sn.sched == SchedStrategy.INF);
T = V .* X;
if isUpper
    R = (1 ./ sn.rates(:)) * N;
    R(isinf_) = 1 ./ sn.rates(isinf_);
    C = Zt + N*sum(D);
else
    R = 1 ./ sn.rates(:);
    C = Zt + sum(D);
end
Q = T .* R;
% Utilization law per SERVER: without the nservers divisor a multiserver
% station reports U > 1 (ssd/ldbcmp/auto reach here)
U = T ./ (max(1, sn.nservers(:)) .* sn.rates(:));
U(isinf_) = Q(isinf_);
lG = - N*log(X);
% Return the scalar bound as XN too: leaving it unset makes every ba_fill
% family report an EMPTY system throughput, and the AUTO composite -- which
% probes candidates by their X -- then discards each of them silently
Xout = X;
end

% Map a SchedStrategy id to a Majumdar-Woodside discipline code:
% 0=FIFO, 1=PS, 2=non-preemptive priority, 3=preemptive priority,
% 4=ABA full-contention (discipline-independent).
function code = mwrbb_disc_code(s)
switch s
    case {SchedStrategy.FCFS}
        code = 0;   % FIFO (Theorem 1 / Lemma 1)
    case {SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS, ...
          SchedStrategy.PSPRIO, SchedStrategy.DPSPRIO, SchedStrategy.GPSPRIO}
        code = 1;   % processor sharing (Lemma 2)
    case {SchedStrategy.HOL}   % HOL == FCFSPRIO: non-preemptive priority
        code = 2;   % non-preemptive priority (Lemmas 4-5)
    case {SchedStrategy.FCFSPRPRIO, SchedStrategy.LCFSPRPRIO}
        code = 3;   % preemptive-resume priority (Lemma 3)
    otherwise
        code = 4;   % non-FCFS work-conserving: ABA full-contention bound
end
end


% Build the MAP-AMVA LP parameters of Casale-Smirni (DSN 2009) from SN.
%
% The LP carries phases at ONE queue and requires it to be the LAST: q(i,j,k,h)
% reads the scalar muM(i) for i < M and the (D0,D1) pair muMAP/v for i == M. A
% model whose phase-carrying station sits elsewhere is therefore PERMUTED, not
% refused; PERM is that permutation and the caller inverts it to report.
% VP and SP are the visit ratios and mean service times in the same permuted
% order, which is what the utilization law is applied in.
function [params, perm, Vp, Sp, N] = ba_mapamva_params(sn)
M = sn.nstations;
N = sn.njobs(1);
PH = sn.proc;

if any(sn.sched == SchedStrategy.INF)
    line_error(mfilename, ['Method ''mapamva'' does not support delay (infinite-server) ' ...
        'stations: the MAP-AMVA program of Casale-Smirni (DSN 2009) is written for a ' ...
        'network of queues and the paper names the delay extension as open work. Use a ' ...
        'QRF method, which carries the load-dependent rate law.']);
end

% Phase order per station. One phase is an exponential server, which enters the
% LP as the scalar rate muM(i); more than one is the (D0,D1) pair, which only
% queue M can hold.
Kph = ones(M,1);
for i = 1:M
    if ~isempty(PH{i}{1})
        Kph(i) = size(PH{i}{1}{1}, 1);
    end
end
phased = find(Kph > 1);
if numel(phased) > 1
    line_error(mfilename, ['Method ''mapamva'' carries phases at ONE station: the LP gives ' ...
        'queue M the (D0,D1) pair and every other queue a scalar rate. Stations %s are all ' ...
        'non-exponential. Use a QRF method, whose q carries a phase at every station.'], ...
        mat2str(phased(:)'));
end
if isempty(phased)
    % Every station exponential: the program is still the right one, it just
    % degenerates to K = 1, where the per-phase variables collapse and the
    % balances become the product-form ones of mapqn_bnd_lr_pf.
    mapIdx = M;
else
    mapIdx = phased;
end
perm = [setdiff(1:M, mapIdx), mapIdx];
K = Kph(mapIdx);

if isempty(PH{mapIdx}{1})
    D0 = -1; D1 = 1;
else
    D0 = PH{mapIdx}{1}{1};
    D1 = PH{mapIdx}{1}{2};
end
% muMAP(k,h) is the completion rate out of phase k landing in phase h, i.e.
% D1(k,h); v(k,h) is the background phase change that completes no job, i.e.
% D0 off the diagonal. Same (from,to) convention as solver_ba_qrf_analyzer --
% writing either as its transpose is invisible for a reversible D0 and silently
% reverses the phase order of an Erlang.
muMAP = D1;
v = D0;
v(1:(K+1):end) = 0;

muM = zeros(max(M-1,0),1);
for a = 1:(M-1)
    muM(a) = sn.rates(perm(a),1);
end

r = zeros(M);
for a = 1:M
    for b = 1:M
        r(a,b) = sn.rt(perm(a), perm(b));
    end
end

V = sn.visits{1}(:);
stimes = zeros(M,1);
for i = 1:M
    if ~isempty(PH{i}{1})
        stimes(i) = map_mean(PH{i}{1});
    end
end
Vp = V(perm);
Sp = stimes(perm);

params = struct();
params.M = M;
params.N = N;
params.K = K;
params.muM = muM;
params.muMAP = muMAP;
params.r = r;
params.v = v;
params.verbose = false;
end
