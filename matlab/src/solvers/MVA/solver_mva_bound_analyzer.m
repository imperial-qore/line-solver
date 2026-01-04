function [Q,U,R,T,C,X,lG,runtime,iter] = solver_mva_bound_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,ITER] = SOLVER_MVA_BOUND_ANALYZER(QN, OPTIONS)

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
end
end
