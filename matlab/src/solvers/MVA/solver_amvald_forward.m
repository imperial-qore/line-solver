function [Wchain, STeff] = solver_amvald_forward(M, K, nservers, schedparam, lldscaling, cdscaling, ljdscaling, ljdcutoffs, ljcdscaling, ljcdcutoffs, sched, classprio, gamma, tau, Qchain_in, Xchain_in, Uchain_in, STchain_in, Vchain_in, Nchain_in, SCVchain_in, options)

%Uhiprio = zeros(M,K); % utilization due to "permanent jobs" in DPS

if isempty(gamma)
    gamma = zeros(M,K);
end

Nt = sum(Nchain_in(isfinite(Nchain_in)));
if all(isinf(Nchain_in))
    delta = 1;
else
    delta = (Nt - 1) / Nt;
end

deltaclass = (Nchain_in - 1) ./ Nchain_in;
deltaclass(isinf(Nchain_in)) = 1;

ocl = find(isinf(Nchain_in));
ccl = find(isfinite(Nchain_in) & Nchain_in>0);
nnzclasses = find(Nchain_in>0);
nnzclasses_eprio = cell(1,length(nnzclasses));
nnzclasses_hprio = cell(1,length(nnzclasses));
nnzclasses_ehprio = cell(1,length(nnzclasses));
if max(classprio) ~= min(classprio)
    for r = nnzclasses
        nnzclasses_eprio{r} = intersect(nnzclasses, find(classprio == classprio(r))); % equal prio
        nnzclasses_hprio{r} = intersect(nnzclasses, find(classprio < classprio(r))); % higher prio (lower value = higher priority)
        nnzclasses_ehprio{r} = intersect(nnzclasses, find(classprio <= classprio(r))); % equal or higher prio
    end
else
    for r = nnzclasses
    nnzclasses_eprio{r} = nnzclasses;
    nnzclasses_hprio{r} = [];
    nnzclasses_ehprio{r} = [];
    end
end

%% evaluate lld and cd correction factors
totArvlQlenSeenByOpen = zeros(K,M,1);
interpTotArvlQlen = zeros(M,1);
totArvlQlenSeenByClosed = zeros(M,K);
stationaryQlen = zeros(M,K);
selfArvlQlenSeenByClosed = zeros(M,K);
for k=1:M
    interpTotArvlQlen(k) = delta * sum(Qchain_in(k,nnzclasses));
    for r = nnzclasses
        selfArvlQlenSeenByClosed(k,r) = deltaclass(r) * Qchain_in(k,r); % qlen of same class as arriving one
        switch sched(k)
            case {SchedStrategy.HOL}
                totArvlQlenSeenByOpen(r,k) = sum(Qchain_in(k,nnzclasses_ehprio{r}));
                totArvlQlenSeenByClosed(k,r) = deltaclass(r) * Qchain_in(k,r) + sum(Qchain_in(k,setdiff(nnzclasses_ehprio{r},r)));
            otherwise
                totArvlQlenSeenByOpen(r,k) = sum(Qchain_in(k,nnzclasses));
                totArvlQlenSeenByClosed(k,r) = deltaclass(r) * Qchain_in(k,r) + sum(Qchain_in(k,nnzclasses)) - Qchain_in(k,r);
        end
        stationaryQlen(k,r) = Qchain_in(k,r); % qlen of same class as arriving one
    end
end

%% all methods use limited-load dependent (LLD) and class-dependent (CD) corrections from QD-AMVA
if isempty(lldscaling)
    lldscaling = ones(M,ceil(Nt));
end

switch options.method
    case {'lin', 'qdlin'}
        if isempty(nnzclasses)
            lldterm = pfqn_lldfun(1 + interpTotArvlQlen, lldscaling);
        else
            lldterm = ones(M,K);
            for r=nnzclasses
                if ~isempty(ccl)
                    lldterm(:,r) = pfqn_lldfun(1 + interpTotArvlQlen + Nchain_in(ccl)*permute(gamma(r,k,ccl),3:-1:1) - gamma(r,k,r), lldscaling);
                end
            end
        end
    otherwise
        lldterm = pfqn_lldfun(1 + interpTotArvlQlen, lldscaling);
end

cdterm = ones(M,K);
switch options.method
    case {'lin', 'qdlin'}
        for r=nnzclasses
            if ~isempty(cdscaling)
                if isfinite(Nchain_in(r))
                    cdterm(:,r) = pfqn_cdfun(1 + selfArvlQlenSeenByClosed + (Nchain_in(r) - 1)*gamma(r,k,r), cdscaling); % qd-amva class-dependence term
                else
                    cdterm(:,r) = pfqn_cdfun(1 + stationaryQlen + (Nchain_in(r) - 1)*gamma(r,k,r), cdscaling); % qd-amva class-dependence term
                end
            end
        end
    otherwise
        for r=nnzclasses
            if ~isempty(cdscaling)
                if isfinite(Nchain_in(r))
                    cdterm(:,r) = pfqn_cdfun(1 + selfArvlQlenSeenByClosed, cdscaling); % qd-amva class-dependence term
                else
                    cdterm(:,r) = pfqn_cdfun(1 + stationaryQlen, cdscaling); % qd-amva class-dependence term
                end
            end
        end
end

%% Joint dependence correction terms
% Both LJD and LJCD are indexed by population vector n = (n1, n2, ..., nK)
% - LJD: scaling = table[idx(n)] - same scaling applied to all classes
% - LJCD: scaling_r = table_r[idx(n)] - each class r has its own table
% If LJCD is set for a station, it overrides LJD for that station.
%
% For FES (Flow-Equivalent Server), the scaling tables contain throughputs X(n).
% Since these are rates, we divide to get time scaling: ljdterm = 1/X(n).
% For multi-class, we interpolate using neighboring population states.

ljdterm = ones(M, K);  % Combined term: (M x K) matrix
hasLjcd = ~isempty(ljcdscaling) && any(~cellfun(@isempty, ljcdscaling));
hasLjd = ~isempty(ljdscaling) && any(~cellfun(@isempty, ljdscaling));


for k = 1:M
    % Get population vector from queue lengths at station k
    nvec = Qchain_in(k, :);

    if hasLjcd && ~isempty(ljcdscaling{k})
        % LJCD: per-class scaling
        % The scaling table contains throughputs X(n) indexed by population n.
        % For FES, the "population" is the number of jobs at the station.
        % Use ceiling of queue length + 1 to account for the arriving job.
        cutoffs = ljcdcutoffs{k};

        % For FES, use ceiling of queue length as the population index
        % The +1 helps for single-class but not multi-class
        if K == 1
            nIndex = ceil(nvec + 1);
            nIndex = max(1, min(nIndex, cutoffs));
        else
            nIndex = ceil(nvec);
            nIndex = max(0, min(nIndex, cutoffs));
        end
        idx = ljd_linearize(nIndex, cutoffs);

        for r = 1:K
            if r <= length(ljcdscaling{k}) && ~isempty(ljcdscaling{k}{r})
                table = ljcdscaling{k}{r};

                % Lookup throughput at population index
                if idx <= length(table)
                    Xval = table(idx);
                else
                    Xval = 0;
                end

                % Convert rate to time scaling (divide by throughput)
                if Xval > GlobalConstants.FineTol
                    ljdterm(k, r) = 1.0 / Xval;
                else
                    ljdterm(k, r) = 1.0;  % No scaling when throughput is zero
                end

            end
        end
    elseif hasLjd && ~isempty(ljdscaling{k})
        % LJD: single scaling applied to all classes
        % Lookup: idx = linearize(n1, ..., nK), then table[idx]
        nClamped = min(nvec, ljdcutoffs(k, :));
        idx = ljd_linearize(nClamped, ljdcutoffs(k, :));
        if idx <= length(ljdscaling{k})
            ljdterm(k, :) = ljdscaling{k}(idx);  % Same value for all classes
        end
    end
end

switch options.config.multiserver
    case 'softmin' % softmin on every station, may not always converge
        switch options.method
            case {'lin', 'qdlin'}
                g = 0;
                for r=ccl
                    g = g + ((Nt-1)/Nt) * Nchain_in(r) * gamma(ccl,:,r);
                end
                msterm = pfqn_lldfun(1 + interpTotArvlQlen + mean(g, 1)', [], nservers); % if native qd then account for multiserver in the correction terms
            otherwise
                msterm = pfqn_lldfun(1 + interpTotArvlQlen + (Nt-1)*mean(gamma(ccl,:), 1)', [], nservers); % if native qd then account for multiserver in the correciton terms
        end
    case 'seidmann'
        msterm = ones(M,1) ./ nservers(:);
        msterm(msterm==0) = 1; % infinite server case
    case 'default' % softmin on ps stations, seidmann on fcfs stations
        switch options.method
            case {'lin', 'qdlin'}
                g = 0;
                for r=ccl
                    g = g + ((Nt-1)/Nt) * Nchain_in(r) * gamma(ccl,:,r);
                end
                msterm = pfqn_lldfun(1 + interpTotArvlQlen + mean(g,1)', [], nservers); % if native qd then account for multiserver in the correction terms
            otherwise
                g = 0;
                for r=ccl
                    g = g + (Nt-1)*gamma(r,:);
                end
                msterm = pfqn_lldfun(1 + interpTotArvlQlen + mean(g), [], nservers); % if native qd then account for multiserver in the correction terms
        end        
        fcfstypeset = find(sched==SchedStrategy.FCFS | sched==SchedStrategy.SIRO | sched==SchedStrategy.LCFSPR);
        msterm(fcfstypeset) = 1 ./ nservers(fcfstypeset);
    otherwise
        line_error(mfilename,'Unrecognized multiserver approximation method');
end

Wchain = zeros(M,K);

STeff = zeros(size(STchain_in)); % effective service time
lldterm = repmat(lldterm, 1, K);
% ljdterm is already (M, K) - contains either LJD or LJCD values
for r=nnzclasses
    for k=1:M
        STeff(k,r) = STchain_in(k,r) * lldterm(k,r) * msterm(k) * cdterm(k,r) * ljdterm(k,r);
    end
end

%% if amva.qli or amva.fli, update now totArvlQlenSeenByClosed with STeff
switch options.method
    case 'qli' % Wang-Sevcik queue line
        infset = sched == SchedStrategy.INF;
        for k=1:M
            switch sched(k)
                case {SchedStrategy.HOL}
                    for r = nnzclasses
                        if Nchain_in(r) == 1
                            totArvlQlenSeenByClosed(k,r) = sum(Qchain_in(k,nnzclasses_ehprio{r})) - Qchain_in(k,r);
                        else
                            qlinum = STeff(k,r) * (1+sum(Qchain_in(k,nnzclasses_ehprio{r})) - Qchain_in(k,r));
                            qliden = sum(STeff(infset,r));
                            for m=1:M
                                qliden = qliden + STeff(m,r) * (1+sum(Qchain_in(m,nnzclasses_ehprio{r})) - Qchain_in(m,r));
                            end
                            totArvlQlenSeenByClosed(k,r) = sum(Qchain_in(k,nnzclasses_ehprio{r})) - (1/(Nchain_in(r)-1))*(Qchain_in(k,r) - qlinum/qliden);
                        end
                    end
                otherwise
                    for r = nnzclasses
                        if Nchain_in(r) == 1
                            totArvlQlenSeenByClosed(k,r) = sum(Qchain_in(k,nnzclasses)) - Qchain_in(k,r);
                        else
                            qlinum = STeff(k,r) * (1+sum(Qchain_in(k,nnzclasses)) - Qchain_in(k,r));
                            qliden = sum(STeff(infset,r));
                            for m=1:M
                                qliden = qliden + STeff(m,r) * (1+sum(Qchain_in(m,nnzclasses)) - Qchain_in(m,r));
                            end
                            totArvlQlenSeenByClosed(k,r) = sum(Qchain_in(k,nnzclasses)) - (1/(Nchain_in(r)-1))*(Qchain_in(k,r) - qlinum/qliden);
                        end
                    end
            end
        end
    case 'fli' % Wang-Sevcik fraction line
        infset = sched == SchedStrategy.INF;
        for k=1:M
            switch sched(k)
                case {SchedStrategy.HOL}
                    for r = nnzclasses
                        if Nchain_in(r) == 1
                            totArvlQlenSeenByClosed(k,r) = sum(Qchain_in(k,nnzclasses_ehprio{r})) - Qchain_in(k,r);
                        else
                            qlinum = STeff(k,r) * (1+sum(Qchain_in(k,nnzclasses_ehprio{r})) - Qchain_in(k,r));
                            qliden = sum(STeff(infset,r));
                            for m=1:M
                                qliden = qliden + STeff(m,r) * (1+sum(Qchain_in(m,nnzclasses_ehprio{r})) - Qchain_in(m,r));
                            end
                            totArvlQlenSeenByClosed(k,r) = sum(Qchain_in(k,nnzclasses_ehprio{r})) - (2/Nchain_in(r))*Qchain_in(k,r) + qlinum/qliden;
                        end
                    end
                otherwise
                    for r = nnzclasses
                        if Nchain_in(r) == 1
                            totArvlQlenSeenByClosed(k,r) = sum(Qchain_in(k,nnzclasses)) - Qchain_in(k,r);
                        else
                            qlinum = STeff(k,r) * (1+sum(Qchain_in(k,nnzclasses)) - Qchain_in(k,r));
                            qliden = sum(STeff(infset,r));
                            for m=1:M
                                qliden = qliden + STeff(m,r) * (1+sum(Qchain_in(m,nnzclasses)) - Qchain_in(m,r));
                            end
                            totArvlQlenSeenByClosed(k,r) = sum(Qchain_in(k,nnzclasses)) - (2/Nchain_in(r))*Qchain_in(k,r) + qlinum/qliden;
                        end
                    end

            end
        end
end

%% compute response times from current queue-lengths
for ir=1:length(nnzclasses)
    r=nnzclasses(ir);
    sd = nnzclasses;
    sd(ir)=[];
    %sdprio = setdiff(nnzclasses_ehprio{r},r);

    for k=1:M

        switch sched(k)
            case SchedStrategy.INF
                Wchain(k,r) = STeff(k,r);

            case SchedStrategy.PS
                switch options.method
                    case {'qd', 'lin', 'qdlin'} % QD-LIN interpolation
                        switch options.config.multiserver
                            case 'seidmann' % in this case, qd handles only lld and cd scalings
                                Wchain(k,r) = STeff(k,r) * (nservers(k)-1); % multi-server correction with serial think time, (1/nservers(k)) term already in STeff
                                if ismember(r,ocl)
                                    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * (1 + totArvlQlenSeenByOpen(r,k));
                                else
                                    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * (1 + interpTotArvlQlen(k) + Nchain_in(ccl)*permute(gamma(r,k,ccl),3:-1:1) - gamma(r,k,r));
                                end
                            case {'default','softmin'}
                                if ismember(r,ocl)
                                    Wchain(k,r) = STeff(k,r) * (1 + totArvlQlenSeenByOpen(r,k));
                                else
                                    switch options.method
                                        case {'lin', 'qdlin'} % Linearizer
                                            Wchain(k,r) = Wchain(k,r) + STeff(k,r) * (1 + interpTotArvlQlen(k) + Nchain_in(ccl)*permute(gamma(r,k,ccl),3:-1:1) - gamma(r,k,r));
                                        otherwise
                                            Wchain(k,r) = STeff(k,r) * (1 + interpTotArvlQlen(k) + (Nt-1)*gamma(r,k));
                                    end
                                end
                        end
                    otherwise 
                        switch options.config.multiserver
                            case 'seidmann'
                                Wchain(k,r) = STeff(k,r) * (nservers(k)-1); % multi-server correction with serial think time, (1/nservers(k)) term already in STeff
                                if ismember(r,ocl)
                                    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * (1 + totArvlQlenSeenByOpen(r,k));
                                else
                                    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * (1 + totArvlQlenSeenByClosed(k) + (Nt-1)*gamma(r,k));
                                end
                            case {'default','softmin'}
                                if ismember(r,ocl)
                                    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * (1 + totArvlQlenSeenByOpen(r,k));
                                else
                                    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * (1 + totArvlQlenSeenByClosed(k) + (Nt-1)*gamma(r,k));
                                end
                        end
                end

            case SchedStrategy.DPS
                if nservers(k)>1
                    line_error(mfilename,'Multi-server DPS not supported yet in AMVA solver.')
                else
                    w = schedparam; % DPS weight
                    tss = Inf; % time-scale separation threshold, this was originally at 5, presently it is disabled
                    %Uhiprio(k,r) = sum(Uchain(k,w(k,:)>tss*w(k,r))); % disabled at the moment
                    %STeff(k,r) = STeff(k,r) / (1-Uhiprio(k,r));

                    Wchain(k,r) = STeff(k,r) * (1 + selfArvlQlenSeenByClosed(k,r)); % class-r
                    for s=sd  % slowdown due to classes s!=r
                        if w(k,s) == w(k,r) % handle gracefully 0/0 case
                            Wchain(k,r) = Wchain(k,r) + STeff(k,r) * stationaryQlen(k,s);
                        elseif w(k,s)/w(k,r)<=tss
                            Wchain(k,r) = Wchain(k,r) + STeff(k,r) * stationaryQlen(k,s) * w(k,s)/w(k,r);
                        elseif w(k,s)/w(k,r)>tss
                            % if there is time-scale separation, do nothing
                            % all is accounted for by 1/(1-Uhiprio)
                        end
                    end
                end

            case {SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.LCFSPR}
                if STeff(k,r) > 0
                    Uchain_r = Uchain_in ./ repmat(Xchain_in,M,1) .* (repmat(Xchain_in,M,1) + repmat(tau(r,:),M,1));

                    if nservers(k)>1
                        deltaclass_r = ones(size(Xchain_in));
                        deltaclass_r(r) = deltaclass(r);
                        if sum(deltaclass_r .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:)) < 0.75 % light-load case
                            Bk = ((deltaclass_r .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:))); % note: this is in 0-1 as a utilization
                        else % high-load case
                            Bk = (deltaclass_r .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:)).^(nservers(k)-1);
                        end
                    else
                        Bk = ones(1,K);
                    end

                    hasLjd = ~isempty(ljdscaling) && any(~cellfun(@isempty, ljdscaling));
                    if nservers(k)==1 && (~isempty(lldscaling) || ~isempty(cdscaling) || hasLjd)
                        switch options.config.highvar % high SCV
                            case 'hvmva'
                                Wchain(k,r) = STeff(k,r) * (1-sum(Uchain_r(k,ccl)));
                                for s=ccl
                                    Wchain(k,r) = Wchain(k,r) + STeff(k,s) * Uchain_r(k,s) * (1 + SCVchain_in(k,s))/2; % high SCV
                                end
                            otherwise % default
                                Wchain(k,r) = STeff(k,r);
                        end
                        if any(ismember(ocl,r))
                            Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * stationaryQlen(k,r) + STeff(k,sd)*stationaryQlen(k,sd)');
                        else
                            switch options.method
                                %case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                %    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r) + STeff(k,sd)*stationaryQlen(k,sd)') + (STeff(k,ccl).*Nchain(ccl)*permute(gamma(r,k,ccl),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                %case {'default', 'amva.lin', 'lin', 'amva.qdlin', 'qdlin'}
                                %    fraction = Qchain_in ./ repmat(Nchain_in, size(Qchain_in, 1), 1);
                                %    Wchain(k, r) = Wchain(k, r) + STeff(k, r) * (Nchain_in(ccl) * fraction(k, ccl)' - fraction(k, r)) + STeff(k, r) * (Nchain_in(ccl)*permute(gamma(r,k,ccl),3:-1:1) - gamma(r,k,r));
                                otherwise
                                    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r) + STeff(k,sd)*stationaryQlen(k,sd)');
                            end
                        end
                    else
                        switch options.config.multiserver
                            case 'softmin'
                                Wchain(k,r) = STeff(k,r); % high SCV
                                if any(ismember(ocl,r))
                                    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * stationaryQlen(k,r) * Bk(r) + STeff(k,sd) * (stationaryQlen(k,sd) .* Bk(sd))';
                                else
                                    % FCFS approximation + reducing backlog proportionally to server utilizations; somewhat similar to
                                    % Rolia-Sevcik -  method of layers - Sec IV.
                                    switch options.method
                                        %case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                        %    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) + STeff(k,sd) * (stationaryQlen(k,sd) .* Bk(sd))' + (STeff(k,ccl).*Nchain(ccl)*permute(gamma(r,k,ccl),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                        otherwise
                                            Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) + STeff(k,sd) * (stationaryQlen(k,sd) .* Bk(sd))';
                                    end
                                end
                            case {'default','seidmann'}
                                Wchain(k,r) = STeff(k,r) * (nservers(k)-1); % multi-server correction with serial think time, (1/nservers(k)) term already in STeff
                                Wchain(k,r) = Wchain(k,r) + STeff(k,r); % high SCV
                                if any(ismember(ocl,r))
                                    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * deltaclass(r) * stationaryQlen(k,r)*Bk(r) + STeff(k,sd).*Bk(sd)*stationaryQlen(k,sd)');
                                else
                                    % FCFS approximation + reducing backlog proportionally to server utilizations; somewhat similar to
                                    % Rolia-Sevcik -  method of layers - Sec IV. (1/nservers(k)) term already in STeff
                                    switch options.method
                                        %case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                        %    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r) + STeff(k,sd).*Bk(sd)*stationaryQlen(k,sd)') + (STeff(k,ccl).*Nchain(ccl)*permute(gamma(r,k,ccl),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                        otherwise
                                            Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r) + STeff(k,sd).*Bk(sd)*stationaryQlen(k,sd)');
                                    end
                                end
                        end
                    end
                end

            case {SchedStrategy.HOL} % non-preemptive priority
                if STeff(k,r) > 0
                    Uchain_r = Uchain_in ./ repmat(Xchain_in,M,1) .* (repmat(Xchain_in,M,1) + repmat(tau(r,:),M,1));

                    switch options.config.np_priority
                        case {'default','cl'} % Chandy-Lakshmi
                            UHigherPrio=0;
                            for h=nnzclasses_hprio{r}
                                UHigherPrio = UHigherPrio + Vchain_in(k,h)*STeff(k,h)*(Xchain_in(h)-Qchain_in(k,h)*tau(h));
                            end
                            prioScaling = min([max([options.tol,1-UHigherPrio]),1-options.tol]);
                        case 'shadow' % Sevcik's shadow server
                            UHigherPrio=0;
                            for h=nnzclasses_hprio{r}
                                UHigherPrio = UHigherPrio + Vchain_in(k,h)*STeff(k,h)*Xchain_in(h);
                            end
                            prioScaling = min([max([options.tol,1-UHigherPrio]),1-options.tol]);
                    end

                    if nservers(k)>1
                        if sum(deltaclass .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:)) < 0.75 % light-load case
                            switch options.config.multiserver
                                case 'softmin'
                                    Bk = ((deltaclass .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:))); % note: this is in 0-1 as a utilization
                                case {'default','seidmann'}
                                    Bk = ((deltaclass .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:)) / nservers(k)); % note: this is in 0-1 as a utilization
                            end
                        else % high-load case
                            switch options.config.multiserver
                                case 'softmin'
                                    Bk = (deltaclass .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:)).^nservers(k); % Rolia
                                case {'default','seidmann'}
                                    Bk = (deltaclass .* Xchain_in .* Vchain_in(k,:) .* STeff(k,:) / nservers(k)).^nservers(k); % Rolia
                            end
                        end
                    else
                        Bk = ones(1,K);
                    end

                    hasLjd = ~isempty(ljdscaling) && any(~cellfun(@isempty, ljdscaling));
                    if nservers(k)==1 && (~isempty(lldscaling) || ~isempty(cdscaling) || hasLjd)
                        switch options.config.highvar % high SCV
                            case 'hvmva'
                                Wchain(k,r) = (STeff(k,r) / prioScaling) * (1-sum(Uchain_r(k,ccl)));
                                for s=ccl
                                    UHigherPrio_s=0;
                                    for h=nnzclasses_hprio{s}
                                        UHigherPrio_s = UHigherPrio_s + Vchain_in(k,h)*STeff(k,h)*(Xchain_in(h)-Qchain_in(k,h)*tau(h));
                                    end
                                    prioScaling_s = min([max([options.tol,1-UHigherPrio_s]),1-options.tol]);
                                    Wchain(k,r) = Wchain(k,r) + (STeff(k,s) / prioScaling_s) * Uchain_r(k,s) * (1 + SCVchain_in(k,s))/2; % high SCV correction
                                end
                            otherwise % default
                                Wchain(k,r) = STeff(k,r) / prioScaling;
                        end

                        if any(ismember(ocl,r))
                            Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * stationaryQlen(k,r)) / prioScaling;
                        else
                            %switch options.method
                            %case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                            %    %Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r) + STeff(k,sdprio)*stationaryQlen(k,sdprio)') + (STeff(k,[r,sdprio]).*Nchain([r,sdprio])*permute(gamma(r,k,[r,sdprio]),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                            %    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r) - STeff(k,r)*gamma(r,k,r)) / prioScaling;
                            %    otherwise
                            %Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r) + STeff(k,sdprio)*stationaryQlen(k,sdprio)');
                            Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r)) / prioScaling;
                            %end
                        end
                    else
                        switch options.config.multiserver
                            case 'softmin'
                                Wchain(k,r) = STeff(k,r) / prioScaling; % high SCV
                                if any(ismember(ocl,r))
                                    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * stationaryQlen(k,r) * Bk(r) / prioScaling;
                                else
                                    % FCFS approximation + reducing backlog proportionally to server utilizations; somewhat similar to
                                    % Rolia-Sevcik -  method of layers - Sec IV.
                                    %switch options.method
                                    %case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                    %    %Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) + STeff(k,sdprio) * (stationaryQlen(k,sdprio) .* Bk(sdprio))' + (STeff(k,[r,sdprio]).*Nchain([r,sdprio])*permute(gamma(r,k,[r,sdprio]),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                    %    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) / prioScaling + (STeff(k,[r]).*Nchain([r])*permute(gamma(r,k,[r]),3:-1:1) - STeff(k,r)*gamma(r,k,r)) / prioScaling;
                                    %    otherwise
                                    %Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) + STeff(k,sdprio) * (stationaryQlen(k,sdprio) .* Bk(sdprio))';
                                    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r) * Bk(r) / prioScaling;
                                    %end
                                end
                            case {'default','seidmann'}
                                Wchain(k,r) = STeff(k,r) * (nservers(k)-1)/prioScaling; % multi-server correction with serial think time, (1/nservers(k)) term already in STeff
                                Wchain(k,r) = Wchain(k,r) + STeff(k,r)/prioScaling; % high SCV
                                if any(ismember(ocl,r))
                                    %Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * stationaryQlen(k,r)*Bk(r) + STeff(k,sdprio).*Bk(sdprio)*stationaryQlen(k,sdprio)');
                                    Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * stationaryQlen(k,r)*Bk(r))/prioScaling;
                                else
                                    % FCFS approximation + reducing backlog proportionally to server utilizations; somewhat similar to
                                    % Rolia-Sevcik -  method of layers - Sec IV. (1/nservers(k)) term already in STeff
                                    %switch options.method
                                    %case {'default', 'amva.lin', 'lin', 'amva.qdlin','qdlin'} % Linearizer
                                    %    %Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r) + STeff(k,sdprio).*Bk(sdprio)*stationaryQlen(k,sdprio)') + (STeff(k,[r,sdprio]).*Nchain([r,sdprio])*permute(gamma(r,k,[r,sdprio]),3:-1:1) - STeff(k,r)*gamma(r,k,r));
                                    %    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r)/prioScaling + (STeff(k,r).*Nchain(r)*permute(gamma(r,k,r),3:-1:1) - STeff(k,r)*gamma(r,k,r))/prioScaling;
                                    %    otherwise
                                    %Wchain(k,r) = Wchain(k,r) + (STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r) + STeff(k,sdprio).*Bk(sdprio)*stationaryQlen(k,sdprio)');
                                    Wchain(k,r) = Wchain(k,r) + STeff(k,r) * selfArvlQlenSeenByClosed(k,r)*Bk(r)/prioScaling;
                                    %end
                                end
                        end
                    end
                end

        end
    end
end
end
