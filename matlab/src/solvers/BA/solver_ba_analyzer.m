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

switch options.method
    case {'lr.upper','lr.lower'}
        % LP-based Linear Reduction bound via mapqn_bnd_lr_pf, distinct from
        % 'qrf.mmi.linear'. see _kb/06-solver-catalog.md for rationale
        line_debug('Using LP linear-reduction bound (%s)', options.method);
        if sn.nclasses ~= 1 || sn.nclosedjobs <= 0
            line_error(mfilename,'Method lr supports single-class closed networks only.');
        end
        if any(sn.nservers(sn.sched ~= SchedStrategy.INF)>1)
            line_error(mfilename,'Unsupported method for a model with multi-server stations.');
        end
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
    case 'aba.upper'
        line_debug('Using ABA upper bound');
        if sn.nclasses==1 && sn.nclosedjobs >0 % closed single-class queueing network
            if any(sn.nservers(sn.sched ~= SchedStrategy.INF)>1)
                line_error(mfilename,'Unsupported method for a model with multi-server stations.');
            end
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
        end
        runtime=toc(T0);
    case 'aba.lower'
        line_debug('Using ABA lower bound');
        if sn.nclasses==1 && sn.nclosedjobs >0 % closed single-class queueing network
            if any(sn.nservers(sn.sched ~= SchedStrategy.INF)>1)
                line_error(mfilename,'Unsupported method for a model with multi-server stations.');
            end
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
        end
        runtime=toc(T0);
    case 'bjb.upper'
        line_debug('Using BJB upper bound');
        if sn.nclasses==1 && sn.nclosedjobs >0 % closed single-class queueing network
            if any(sn.nservers(sn.sched ~= SchedStrategy.INF)>1)
                line_error(mfilename,'Unsupported method for a model with multi-server stations.');
            end
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
        end
        runtime=toc(T0);
    case 'bjb.lower'
        line_debug('Using BJB lower bound');
        if sn.nclasses==1 && sn.nclosedjobs >0 % closed single-class queueing network
            if any(sn.nservers(sn.sched ~= SchedStrategy.INF)>1)
                line_error(mfilename,'Unsupported method for a model with multi-server stations.');
            end
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
        end
        runtime=toc(T0);
    case 'pb.upper'
        line_debug('Using PB upper bound');
        if sn.nclasses==1 && sn.nclosedjobs >0 % closed single-class queueing network
            if any(sn.nservers(sn.sched ~= SchedStrategy.INF)>1)
                line_error(mfilename,'Unsupported method for a model with multi-server stations.');
            end
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
        end
        runtime=toc(T0);
    case 'pb.lower'
        line_debug('Using PB lower bound');
        if sn.nclasses==1 && sn.nclosedjobs >0 % closed single-class queueing network
            if any(sn.nservers(sn.sched ~= SchedStrategy.INF)>1)
                line_error(mfilename,'Unsupported method for a model with multi-server stations.');
            end
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
        end
        runtime=toc(T0);    
    case 'sb.upper'
        line_debug('Using SB upper bound');
        if sn.nclasses==1 && sn.nclosedjobs >0 % closed single-class queueing network
            if any(sn.nservers(sn.sched ~= SchedStrategy.INF)>1)
                line_error(mfilename,'Unsupported method for a model with multi-server stations.');
            end
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
            C(1,1) = Z+A1+(N-1)*(A1*A2+A3)/(A1^2+A2);
            X(1,1) = min([1/Dmax,N / (Z+A1+(N-1)*(A1*A2+A3)/(A1^2+A2))]);
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
        end
        runtime=toc(T0);        
    case 'sb.lower'
        line_debug('Using SB lower bound');
        if sn.nclasses==1 && sn.nclosedjobs >0 % closed single-class queueing network
            if any(sn.nservers(sn.sched ~= SchedStrategy.INF)>1)
                line_error(mfilename,'Unsupported method for a model with multi-server stations.');
            end
            if any(sn.sched == SchedStrategy.INF)
                line_error(mfilename,'Unsupported method for a model with infinite-server stations.');
            end
            V = sn.visits{1}(:);
            D = V(sn.sched ~= SchedStrategy.INF) ./ sn.rates(sn.sched ~= SchedStrategy.INF);
            Z = sum(V(sn.sched == SchedStrategy.INF) ./ sn.rates(sn.sched == SchedStrategy.INF));
            N = sn.nclosedjobs;
            AN = sum(D.^N);
            A1 = sum(D);
            C(1,1) = Z+A1+(N-1)*(AN/A1)^(1/(N-1));
            X(1,1) = N / (Z+A1+(N-1)*(AN/A1)^(1/(N-1)));
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
        end
        runtime=toc(T0);
    case 'gb.upper'
        line_debug('Using GB upper bound');
        if sn.nclasses==1 && sn.nclosedjobs >0 % closed single-class queueing network
            if any(sn.nservers(sn.sched ~= SchedStrategy.INF)>1)
                line_error(mfilename,'Unsupported method for a model with multi-server stations.');
            end
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
                    Q(i,1) = X(1,1) * R(i,1);
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
        end
        runtime=toc(T0);
    case 'gb.lower'
        line_debug('Using GB lower bound');
        if sn.nclasses==1 && sn.nclosedjobs >0 % closed single-class queueing network
            if any(sn.nservers(sn.sched ~= SchedStrategy.INF)>1)
                line_error(mfilename,'Unsupported method for a model with multi-server stations.');
            end
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
                    Q(i,1) = X(1,1) * R(i,1);
                else
                    k = k + 1;
                    Q(i,1) = pfqn_qzgblow(D,N,Z,k);
                    R(i,1) = Q(i,1) / XNup / V(i) ;
                end
            end
            U(:,1) = T(:,1) ./ sn.rates;
            U((sn.sched == SchedStrategy.INF),1) = Q((sn.sched == SchedStrategy.INF),1);
            lG = - N*log(X(1,1)); % approx
        end
        runtime=toc(T0);
    case {'mwba.upper','mwba.lower'}
        line_debug(['Using Majumdar-Woodside robust box bound: ' options.method]);
        if sn.nclosedjobs > 0 && ~any(isinf(sn.njobs)) % fully closed multiclass network
            if any(sn.nservers(sn.sched ~= SchedStrategy.INF) > 1)
                line_error(mfilename,'Unsupported method for a model with multi-server stations.');
            end
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
            nq = find(~isdelay);
            Tchain = zeros(sn.nstations,sn.nchains);
            Uchain = zeros(sn.nstations,sn.nchains);
            Qchain = zeros(sn.nstations,sn.nchains);
            for c=1:sn.nchains
                for i=1:sn.nstations
                    Tchain(i,c) = Xchain(c) * Vchain(i,c);
                    Uchain(i,c) = Xchain(c) * Lchain(i,c);   % utilization law
                    if isdelay(i)
                        Qchain(i,c) = Xchain(c) * Lchain(i,c);
                    else
                        kk = find(nq==i);
                        if strcmp(options.method,'mwba.upper')
                            w = Sq(kk,c);                    % no-contention residence
                        else
                            w = Wlo(kk,c);                   % Theorem 1 residence
                        end
                        Qchain(i,c) = Xchain(c) * Vchain(i,c) * w;
                    end
                end
            end
            [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, Qchain, Uchain, [], Tchain, [], Xchain);
            lG = NaN;
        end
        runtime=toc(T0);
    case {'pbh.upper','pbh.lower'}
        line_debug(['Using PBH (Eager-Sevcik) bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn, options.method);
        ba_reject_multiserver(sn, options.method);
        [Xlo,Xhi] = pfqn_pbh(D, N, Zt, ba_level(options));
        up = strcmp(options.method,'pbh.upper');
        [Q,U,R,T,C,lG] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'cbh.upper','cbh.lower'}
        line_debug(['Using CBH (Dowdy) bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn, options.method);
        ba_reject_multiserver(sn, options.method);
        [Xlo,Xhi] = pfqn_cbh(D, N, Zt, ba_level(options));
        up = strcmp(options.method,'cbh.upper');
        [Q,U,R,T,C,lG] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'pbk.upper','pbk.lower'}
        line_debug(['Using PB(k) bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn, options.method);
        ba_reject_multiserver(sn, options.method);
        [Xlo,Xhi] = pfqn_pbk(D, N, Zt, ba_level(options));
        up = strcmp(options.method,'pbk.upper');
        [Q,U,R,T,C,lG] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'bjbk.upper','bjbk.lower'}
        line_debug(['Using BJB(k) bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn, options.method);
        ba_reject_multiserver(sn, options.method);
        [Xlo,Xhi] = pfqn_bjbk(D, N, Zt, ba_level(options));
        up = strcmp(options.method,'bjbk.upper');
        [Q,U,R,T,C,lG] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'ssd.upper','ssd.lower'}
        line_debug(['Using SSD (Suri-Dallery) multiserver bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn, options.method);
        cvec = sn.nservers(sn.sched ~= SchedStrategy.INF);
        [Xlo,Xhi] = pfqn_ssd(D, N, Zt, cvec(:));
        up = strcmp(options.method,'ssd.upper');
        [Q,U,R,T,C,lG] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'cub.upper','mbjb.lower'}
        % cub.upper: Kerola composite upper bound (upper-only in the paper).
        % mbjb.lower: the multiclass Balanced Job Bounds lower bound (Kerola
        % eq. 10) that seeds it; the natural multiclass counterpart of bjb.
        line_debug(['Using multiclass composite/BJB bound: ' options.method]);
        if sn.nclosedjobs <= 0 || any(isinf(sn.njobs))
            line_error(mfilename,'Method ''%s'' supports fully closed networks only.', options.method);
        end
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
            Qchain(:,c) = Uchain(:,c);               % coarse (bound on X is primary)
        end
        [Q,U,R,T,C,X] = sn_deaggregate_chain_results(sn, Lchain, [], STchain, Vchain, alpha, Qchain, Uchain, [], Tchain, [], Xchain);
        lG = NaN;
        runtime=toc(T0);
    case {'sib.upper','sib.lower'}
        line_debug(['Using SIB (Srinivasan) bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn, options.method);
        ba_reject_multiserver(sn, options.method);
        if Zt > 0
            line_error(mfilename,'Method ''%s'' supports Z=0 (no delay station) only; delay needs the SIB Section-3.2 extension.', options.method);
        end
        [Xlo,Xhi] = pfqn_sib(D, N, 0, ba_level(options));
        up = strcmp(options.method,'sib.upper');
        [Q,U,R,T,C,lG] = ba_fill(sn, V, N, up*Xhi+(~up)*Xlo, Zt, D, up);
        runtime=toc(T0);
    case {'ldbcmp.lower'}
        line_debug(['Using LD-BCMP (Anselmi-Cremonesi) lower throughput bound: ' options.method]);
        [V,Zt,D,N] = ba_sc_demands(sn, options.method);
        % Fixed-rate / delay parameterization (Heffes c=0); LD stations are
        % treated as fixed-rate at their limiting demand.
        [Xlo,~,Qhat] = pfqn_ldbcmp(D, N, Zt, zeros(numel(D),1));
        if isnan(Xlo)
            line_error(mfilename,'Method ''%s'' requires the asymptotic regime N >= Qhat (Qhat=%.4f > N=%d).', options.method, Qhat, N);
        end
        [Q,U,R,T,C,lG] = ba_fill(sn, V, N, Xlo, Zt, D, false);
        runtime=toc(T0);
end
end

% Single-class demand extraction for the hierarchical bound cases: returns
% visit vector V, aggregate think time Zt, per-queue demand vector D and the
% closed population N. Rejects non-single-class models with a clear error.
function [V,Zt,D,N] = ba_sc_demands(sn, method)
if sn.nclasses ~= 1 || sn.nclosedjobs <= 0
    line_error('solver_ba_analyzer', ...
        'Method ''%s'' supports single-class closed networks only.', method);
end
V = sn.visits{1}(:);
isinf_ = (sn.sched == SchedStrategy.INF);
Zt = sum(V(isinf_) ./ sn.rates(isinf_));
D  = V(~isinf_) ./ sn.rates(~isinf_);
N  = sn.nclosedjobs;
end

function ba_reject_multiserver(sn, method)
if any(sn.nservers(sn.sched ~= SchedStrategy.INF) > 1)
    line_error('solver_ba_analyzer', ...
        'Method ''%s'' does not support multi-server stations (use ''ssd'').', method);
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
function [Q,U,R,T,C,lG] = ba_fill(sn, V, N, X, Zt, D, isUpper)
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
U = T ./ sn.rates(:);
U(isinf_) = Q(isinf_);
lG = - N*log(X);
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
