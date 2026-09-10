function [RD] = solver_mam_passage_time(sn, PH, options)
% [RD] = SOLVER_MAM_PASSAGE_TIME(QN, PH, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% generate local state spaces
M = sn.nstations;
K = sn.nclasses;
N = sn.njobs';
rt = sn.rt;
V = sn.visits;

% Get number of CDF points from options
if isfield(options, 'config') && isfield(options.config, 'num_cdf_pts')
    n_cdf_pts = options.config.num_cdf_pts;
else
    n_cdf_pts = 100;
end

Q = zeros(M,K);
U = zeros(M,K);
R = zeros(M,K);
T = zeros(M,K);
C = zeros(1,K);
X = zeros(1,K);

if M==2 && all(isinf(N))
    % open queueing system (one node is the external world)
    pie = {};
    S = {};
    for ist=1:M
        switch sn.sched(ist)
            case SchedStrategy.EXT
                na = cellfun(@(x) length(x{1}),{PH{ist}{:}});
                A = {PH{ist}{1}{1},PH{ist}{1}{2},PH{ist}{1}{2}};
                for k=2:K
                    A = mmap_super(A,{PH{ist}{k}{1},PH{ist}{k}{2},PH{ist}{k}{2}});
                end
                idx_arv = ist;
            case {SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO}
                row = size(S,1) + 1;
                for k=1:K
                    PH{ist}{k} = map_scale(PH{ist}{k}, map_mean(PH{ist}{k})/sn.nservers(ist));
                    pie{k} = map_pie(PH{ist}{k});
                    S{k} = PH{ist}{k}{1};
                end
                idx_q = ist;
                is_ps = false;
            case SchedStrategy.PS
                % Processor-sharing queue
                row = size(S,1) + 1;
                for k=1:K
                    % For PS, service times remain exponential (no scaling needed for distribution)
                    S{k} = PH{ist}{k}{1};  % Service process subgenerator
                    pie{k} = map_pie(PH{ist}{k});
                end
                idx_q = ist;
                is_ps = true;
            otherwise
                line_error(mfilename,'Unsupported scheduling strategy');
        end
    end

    % Priorities select the law only under a priority DISCIPLINE: a plain FCFS
    % or PS queue serves in arrival or processor order whatever the classprio
    % column says, exactly as solver_mam_basic gates the same analyzers
    if (sn.sched(idx_q)==SchedStrategy.HOL || sn.sched(idx_q)==SchedStrategy.FCFSPRPRIO) ...
            && any(sn.classprio ~= sn.classprio(1)) % if priorities are not identical
        [uK,iK] = unique(sn.classprio);
        if length(uK) == length(sn.classprio) % if all priorities are different
            % BUTools convention: D1=lowest priority, DK=highest priority
            % LINE convention: lower value = higher priority
            % unique() returns ascending order, so reverse for BUTools
            iK = flipud(iK(:));
            % Neither vendored analyzer exports 'stDistrPH', so the law is
            % TABULATED: one solve for two sojourn moments sizes a shared grid
            % (mean + 10 sigma over the classes), a second tabulates the CDF on
            % it. The [outputs{iK'}] assignment lands each BUTools slot back on
            % its LINE class.
            Moms = cell(1,K);
            if sn.sched(idx_q)==SchedStrategy.FCFSPRPRIO
                [Moms{iK'}] = MMAPPH1PRPR({A{[1;2+iK]}}, {pie{iK}}, {S{iK}}, 'stMoms', 2);
            else % HOL (non-preemptive)
                [Moms{iK'}] = MMAPPH1NPPR({A{[1;2+iK]}}, {pie{iK}}, {S{iK}}, 'stMoms', 2);
            end
            xmax = 0;
            for k=1:K
                m1 = Moms{k}(1);
                sig = sqrt(max(Moms{k}(2) - m1^2, 0));
                xmax = max(xmax, m1 + 10*sig);
            end
            X = linspace(0, xmax, n_cdf_pts);
            % The Erlangization inside stDistr divides by t, so t = 0 is not
            % evaluable; a sojourn time is strictly positive, so F(0) = 0
            % exactly and the grid's first point is prepended, not asked for
            Dist = cell(1,K);
            if sn.sched(idx_q)==SchedStrategy.FCFSPRPRIO
                [Dist{iK'}] = MMAPPH1PRPR({A{[1;2+iK]}}, {pie{iK}}, {S{iK}}, 'stDistr', X(2:end));
            else
                [Dist{iK'}] = MMAPPH1NPPR({A{[1;2+iK]}}, {pie{iK}}, {S{iK}}, 'stDistr', X(2:end));
            end
            RD = cell(M,K);
            for k=1:K
                RD{idx_arv,k} = [];
                RD{idx_q,k} = [[0; Dist{k}(:)], X(:)];
            end
        else
            line_error(mfilename,'SolverMAM requires either identical priorities or all distinct priorities');
        end
    else
        RD = cell(M,K);

        if is_ps
            % MAP/M/1-PS queue: use specialized sojourn time distribution algorithm
            % Check that service processes are exponential (M)
            for k=1:K
                if size(S{k},1) ~= 1
                    line_error(mfilename,'PS queue requires exponential (Markovian) service times');
                end
            end

            % For single-class case, directly compute sojourn time CDF
            if K == 1
                % Extract MAP arrival parameters
                C_map = A{1};  % MAP C matrix
                D_map = A{2};  % MAP D matrix

                % Extract exponential service rate
                mu = -S{1}(1,1);  % Service rate (S is negative generator)

                % Estimate CDF range
                % First get approximate mean from steady-state
                lambda = map_lambda({C_map, D_map});
                rho = lambda / mu;
                approx_mean = 1 / (mu * (1 - rho));  % Approximate M/M/1-PS mean

                % Generate x points (use similar range as FCFS case)
                x_max = approx_mean * 10;  % Conservative upper bound
                x_vals = linspace(0, x_max, n_cdf_pts);

                % Call MAP/M/1-PS sojourn time CDF function
                W_bar = map_m1ps_cdfrespt(C_map, D_map, mu, x_vals);

                % Convert to CDF format (currently W_bar is complementary CDF)
                F = 1 - W_bar(:);

                % Store results
                RD{idx_arv,1} = [];
                RD{idx_q,1} = [F, x_vals(:)];
            else
                % Multi-class PS: aggregate arrivals into single MAP
                % Extract individual service rates
                mu_vec = zeros(1,K);
                for k=1:K
                    mu_vec(k) = -S{k}(1,1);
                end

                % Check if all service rates are equal (required for current implementation)
                if any(abs(mu_vec - mu_vec(1)) > GlobalConstants.FineTol)
                    line_error(mfilename,'Multi-class PS currently requires identical service rates');
                end
                mu = mu_vec(1);

                % Compute aggregated MAP parameters
                C_map = A{1};
                % Sum the CLASS arrival matrices, A{3:end}: A{2} is already the
                % aggregate the classes sum to, so starting the sum there
                % counted every arrival twice and doubled the load
                D_map_sum = sum(cat(3, A{3:end}), 3);

                % Estimate CDF range using aggregated arrival rate
                lambda = map_lambda({C_map, D_map_sum});
                rho = lambda / mu;
                approx_mean = 1 / (mu * (1 - rho));

                x_max = approx_mean * 10;
                x_vals = linspace(0, x_max, n_cdf_pts);

                % Compute sojourn time for aggregated system
                W_bar = map_m1ps_cdfrespt(C_map, D_map_sum, mu, x_vals);
                F = 1 - W_bar(:);

                % Assign same CDF to all classes (approximation)
                for k=1:K
                    RD{idx_arv,k} = [];
                    RD{idx_q,k} = [F, x_vals(:)];
                end
            end
        else
            % FCFS/HOL queue: use existing MMAPPH1FCFS method
            [Ret{1:2*K}] = MMAPPH1FCFS({A{[1,3:end]}}, {pie{:}}, {S{:}}, 'stDistrPH');

            for k=1:K
                alpha{k} = Ret{(k-1)*2+1};
                D0{k} = Ret{(k-1)*2+2};
                RDph{k}={D0{k},(-D0{k})*ones(length(alpha{k}),1)*alpha{k}(:)'};
                % now estimate the range of the CDF
                sigma(k) = sqrt(map_var(RDph{k}));
                mean(k) = map_mean(RDph{k});
                n = 5;
                while map_cdf(RDph{k},mean(k)+n*sigma(k)) < 1-GlobalConstants.FineTol
                    n = n+1;
                end
                % generate CDF points
                X = linspace(0,(mean(k)+n*sigma(k)),n_cdf_pts);
                F = map_cdf(RDph{k},X);
                RD{idx_arv,k} = [];
                RD{idx_q,k} = [F(:),X(:)];
            end
        end
    end
else
    line_warning(mfilename,'This model is not supported by SolverMAM yet. Returning with no result.\n');
    RD = {}; % warn-and-empty, as the JAR does; an unassigned RD errored at the caller
end

end
