function RTret = solver_fluid_passage_time(sn, options)
% RTRET = SOLVER_FLUID_PASSAGE_TIME(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

iter_max = options.iter_max;
tol = options.tol;
y0 = options.init_sol;
stiff = options.stiff;
T0 = 0;
M = sn.nstations;    %number of stations
K = sn.nclasses;    %number of classes
N = sn.nclosedjobs;    %population
Lambda = sn.mu;
Pi = sn.phi;
PH = sn.proc;
rt = sn.rt;
S = sn.nservers;

for j = 1:M
    %Set number of servers in delay station = population
    if isinf(S(j))
        S(j) = N;
    end
end

%% initialization
phases = sn.phases;
slowrate = zeros(M,K);
for j = 1:M
    for k = 1:K
        slowrate(j,k) = Inf;
        slowrate(j,k) = min(slowrate(j,k),min(Lambda{j}{k}(:))); %service completion (exit) rates in each phase
    end
end

%% response time analysis - starting from fixed point found
%stiff = 1;
chains = sn.chains;
nChains = size(chains,1);

RT = [];
for i = 1:sn.nstations
    if sn.nodetype(sn.stationToNode(i)) ~= NodeType.Source
        for k = 1:nChains %once for each chain
            idxClassesInChain = find(chains(k,:)==1);
            for c = idxClassesInChain
                if phases(i,c) > 0
                    [Kc, ode_h_c, y0_c, phases_c, fluid_c] = generate_odes_passagetime(k,c);

                    % determine max integration time
                    nonZeroRates = slowrate(:);
                    nonZeroRates = nonZeroRates( nonZeroRates > tol );
                    T = abs(100/min(nonZeroRates)); % solve ode until T = 100 events with slowest exit rate

                    % indices of new classes at station i
                    idxN = [];
                    for j = i % this used to be at all stations
                        idxN = [idxN sum(sum(phases_c(1:j-1,: ) )) + sum(phases_c(j,1:K)) + [1:sum(phases_c(j,K+1:Kc))] ];
                    end

                    %% ODE analysis
                    fullt = [];
                    fully = [];
                    iter = 1;
                    finished = 0;
                    tref = 0;
                    odeopt = odeset('AbsTol', tol, 'RelTol', tol, 'NonNegative', 1:length(y0_c));
                    while iter <= iter_max && finished == 0
                        trange = [T0, T];
                        % solve ode - ymean_t_iter is the transient solution in stage e
                        try
                            if stiff
                                [t_iter, ymean_t_iter] = ode_solve_stiff(ode_h_c, trange, y0_c, odeopt, options);
                            else
                                [t_iter, ymean_t_iter] = ode_solve(ode_h_c, trange, y0_c, odeopt, options);
                            end
                        catch ME
                            line_printf('\nODE solver failed. Fluid solver switching to default initialization.');
                            odeopt = odeset('AbsTol', tol, 'RelTol', tol, 'NonNegative', 1:length(y0_c));
                            %try
                            [t_iter, ymean_t_iter] = ode_solve(ode_h_c, trange, y0_c, odeopt, options);
                            %catch
                            %    %keyboard
                            %end
                        end
                        %%%
                        iter = iter + 1;
                        if isempty(fullt)
                            fullt = [fullt; t_iter+tref];
                            fully = [fully; ymean_t_iter];
                        else
                            fullt = [fullt; t_iter(2:end)+tref];
                            fully = [fully; ymean_t_iter(2:end,:)];
                        end
                        if sum(ymean_t_iter(end,idxN )) < 10e-10
                            finished = 1;
                        end
                        tref = tref + t_iter(end);
                        y0_c = ymean_t_iter(end,:);
                    end

                    % retrieve response time CDF for class c
                    RT{i,c,1} = fullt;
                    if fluid_c > 0
                        RT{i,c,2} = 1 - sum(fully(:,idxN ),2)/fluid_c;
                    else
                        RT{i,c,2} = ones(size(fullt));
                    end

                    %% Adaptive CDF Refinement - detect and refine large CDF jumps
                    if fluid_c > 0
                        maxCdfJump = 0.001; % Maximum allowed CDF jump (0.1%)
                        cdfValues = RT{i,c,2};
                        for row = 2:length(cdfValues)
                            cdfJump = cdfValues(row) - cdfValues(row-1);
                            if cdfJump > maxCdfJump
                                % Refine this interval
                                t1 = fullt(row-1);
                                t2 = fullt(row);
                                numRefinedPoints = 20;
                                refinedT = linspace(t1, t2, numRefinedPoints)';
                                refinedStates = zeros(numRefinedPoints, size(fully,2));

                                % Integrate to refined time points
                                startState = fully(row-1,:);
                                for rp = 1:numRefinedPoints
                                    try
                                        if stiff
                                            [~, tempState] = ode_solve_stiff(ode_h_c, [t1, refinedT(rp)], startState, odeopt, options);
                                        else
                                            [~, tempState] = ode_solve(ode_h_c, [t1, refinedT(rp)], startState, odeopt, options);
                                        end
                                        refinedStates(rp,:) = max(0, tempState(end,:));
                                    catch
                                        % Linear interpolation fallback
                                        alpha = (refinedT(rp) - t1) / (t2 - t1);
                                        refinedStates(rp,:) = fully(row-1,:) + alpha * (fully(row,:) - fully(row-1,:));
                                    end
                                end

                                % Merge refined points using interpolation
                                newFullt = unique([fullt(1:row-1); refinedT(2:end); fullt(row+1:end)]);
                                newFully = interp1(fullt, fully, newFullt, 'linear', 'extrap');
                                fullt = newFullt;
                                fully = newFully;

                                % Recompute CDF
                                RT{i,c,1} = fullt;
                                RT{i,c,2} = 1 - sum(fully(:,idxN),2)/fluid_c;
                                cdfValues = RT{i,c,2};

                                line_printf('INFO: Added %d refined points between t=%.6f and t=%.6f\n', numRefinedPoints, t1, t2);
                                break; % Only refine first large jump per iteration
                            end
                        end

                        %% Extended Time Interval Logic - extend if first CDF > 1%
                        maxExtendIterations = 10;
                        extendIter = 0;
                        while RT{i,c,2}(1) > 0.01 && extendIter < maxExtendIterations
                            extendIter = extendIter + 1;
                            extendedT = T * (1 + extendIter);

                            % Re-run with extended time
                            try
                                if stiff
                                    [t_ext, y_ext] = ode_solve_stiff(ode_h_c, [T0, extendedT], y0_c, odeopt, options);
                                else
                                    [t_ext, y_ext] = ode_solve(ode_h_c, [T0, extendedT], y0_c, odeopt, options);
                                end
                                fullt = t_ext;
                                fully = y_ext;
                                RT{i,c,1} = fullt;
                                RT{i,c,2} = 1 - sum(fully(:,idxN),2)/fluid_c;
                            catch
                                break; % Stop extending on error
                            end
                        end
                    end

                    if iter > iter_max
                        line_printf('\n');
                        line_warning(mfilename,'Maximum number of iterations reached when computing the response time distribution at station %d in class %d.\n',i,c);
                        line_warning(mfilename,'Response time distributions may be inaccurate. Try increasing option.iter_max (currently at %s).\n',num2str(iter_max));
                    end
                end
            end
        end
    end
end

RTret = {};
for i=1:sn.nstations
    if sn.nodetype(sn.stationToNode(i)) ~= NodeType.Source
        for c=1:sn.nclasses
            RTret{i,c} = [RT{i,c,2},RT{i,c,1}];            
            if ~isempty(RTret{i,c}) && RTret{i,c}(end,1) < 0.995
                line_warning(mfilename,'CDF at station %d in class %d computed only %.3f percent of the total mass.\n',i,c,RTret{i,c}(end,1)*100);
            end
        end
    end
end
return

    function [Kc, ode_h_c, y0_c, phases_c, fluid_c] = generate_odes_passagetime(k,c)
        Kc = K + 1;  % add a single new class
        numTranClasses = Kc - K;
        idxTranCl = zeros(1,K); % indices of the transient class corresponding to each class in the original model for chain k
        idxTranCl(chains(k,:)==1) =  K+1:Kc; % this is just K+1 since we are adding a single new class, but this format may be generalizable
        new_mu = cell(M,1);
        for m=1:M
            new_mu{m,1} = cell(1,Kc);
        end
        new_pi = cell(M,1);
        for m=1:M
            new_pi{m,1} = cell(1,Kc);
        end
        new_rt = zeros(M*Kc, M*Kc); % new routing table
        new_proc = PH;

        for j=1:sn.nstations
            % service rates
            for k=1:K
                new_mu{j}{k} = Lambda{j}{k};
            end
            new_mu{j}{K+1} = Lambda{j}{c};

            % completion probabilities
            for k=1:K
                new_pi{j}{k} = Pi{j}{k};
            end
            new_pi{j}{K+1} = Pi{j}{c};

            % phd distribution
            for r=1:nChains
                new_proc{j}{r} = PH{j}{r};
            end
            new_proc{j}{K+1} = PH{j}{c};
        end

        % routing/switching probabilities
        % among basic classes
        for l = 1:K
            for m = 1:K
                new_rt(l:Kc:end,m:Kc:end) = rt(l:K:end,m:K:end);
            end
        end

        % copy routing table from the original to the transient classes (forward)
        for l = 1:numTranClasses
            for m = 1:numTranClasses
                if sum(sum(rt(c:K:end,idxClassesInChain(m):K:end))) > 0
                    new_rt(K+l:Kc:end,K+m:Kc:end) = rt(c:K:end,idxClassesInChain(m):K:end);
                end
            end
        end

        %phases of transient classes
        phases_c = zeros(M,Kc);
        phases_c(:,1:K) = phases;
        phases_c(:,K+1) = phases(:,c);

        % identify classes in chain that complete
        completingClassesInChain = c;

        %determine final classes (leaves in the class graph)
        for s = completingClassesInChain' % for each completing class
            %routing matrix from a transient class that completes is diverted back into the original classes
            for l = idxClassesInChain
                for j = 1:sn.nstations
                    % return fluid to original class
                    new_rt((i-1)*Kc+idxTranCl(c), (j-1)*Kc+l) = rt((i-1)*K+c, (j-1)*K+l);
                    % delete corresponding transition among transient classes
                    new_rt((i-1)*Kc+idxTranCl(c), (j-1)*Kc+idxTranCl(l)) = 0;
                end
            end
        end

        % setup the ODEs for the new QN
        %        options.method  = 'statedep'; % default doesn't seem to work in some models
        [ode_h_c, ~] = solver_fluid_odes(sn, N, new_mu', new_pi', new_proc, new_rt, S, sn.sched, sn.schedparam, options);

        % setup initial point
        y0_c = zeros(1, sum(sum(phases_c(:,:))));
        fluid_c = 0;
        for j = 1:sn.nstations
            for l = 1:sn.nclasses
                idxNew_jl = sum(sum(phases_c(1:j-1,:))) + sum(phases_c(j,1:l-1));
                idxNew_jt = sum(sum(phases_c(1:j-1,:))) + sum(phases_c(j,1:idxTranCl(l)-1));
                idx_jl = sum(sum(phases(1:j-1,:))) + sum(phases(j,1:l-1));
                if i == j && l==c
                    y0_c( idxNew_jt + 1 ) = sum(y0(idx_jl+1:idx_jl + phases(j,l))); % mass in phases all moved back into phase 1
                    fluid_c = fluid_c + sum(y0(idx_jl+1:idx_jl + phases(j,l)));
                else % leave mass as it is
                    y0_c( idxNew_jl + 1: idxNew_jl + phases_c(j,l)  ) = y0(idx_jl+1:idx_jl + phases(j,l));
                end 
            end
        end
    end

end