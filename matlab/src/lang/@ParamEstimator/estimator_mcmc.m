function estVal = estimator_mcmc(self, nodes)
    % Gibbs Sampling MCMC-based optimization
    % Supports closed, open, and mixed queueing networks.
    % For open classes, uses the open-to-closed equivalence Z_r = N_r / lambda_r.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    % This code is released under the 3-Clause BSD License.

    sn = self.model.getStruct;

    % Get effective population for open classes
    if isfield(self.options, 'openPopulation') && ~isempty(self.options.openPopulation)
        Nopen = self.options.openPopulation;
    else
        Nopen = 100;
    end

    % Determine population per class
    P = zeros(1, sn.nclasses);
    for r = 1:sn.nclasses
        if sn.njobs(r) < Inf
            P(r) = sn.njobs(r);
        else
            P(r) = Nopen;
        end
    end

    % Extract think times
    Z = zeros(1, sn.nclasses);
    allNodes = self.model.getNodes;
    for n = 1:length(allNodes)
        if isa(allNodes{n}, 'Delay')
            svcProc = allNodes{n}.getService;
            for r = 1:sn.nclasses
                if sn.njobs(r) < Inf
                    Z(r) = Z(r) + svcProc{r}.getMean();
                end
            end
        elseif isa(allNodes{n}, 'Source')
            svcProc = allNodes{n}.getService;
            for r = 1:sn.nclasses
                if sn.njobs(r) == Inf
                    lambda_r = 1 / svcProc{r}.getMean();
                    Z(r) = P(r) / lambda_r;
                end
            end
        end
    end

    % obtain per class metrics
    for n=1:size(nodes,2)
        node = nodes{n};

        for r=1:sn.nclasses
            avgAQLen{n,r} = self.getAggrQLen(node); % aggregate queue-length
            if isempty(avgAQLen{n,r})
                error('Transient queue-length data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
            else
                avgAQLen{n,r} = avgAQLen{n,r}.data;
            end
        end
    end

    try
        avgQL = cell2mat(avgAQLen);
        experiments = size(avgQL,1);
        avgQL = reshape(avgQL, size(avgQL,1)/size(avgAQLen, 1), size(avgAQLen, 1), sn.nclasses);
        avgQL = squeeze(mean(avgQL,1))';
    catch me
        switch me.identifier
            case 'MATLAB:catenate:dimensionMismatch'
                error('Sampled metrics have different number of samples, use interpolate() before starting this estimation algorithm.');
        end
    end

    estVal = mcmc_data(avgQL, self.model.getStruct.visits, experiments, self.options.iter_max, P, Z);
end

% mcmc procedure based on the comon data format
function demEst = mcmc_data(avgQL, visits, experiments, ITERMAX, P, Z)

    %% number of resources
    M = size(avgQL,1);
    %% number of classes
    R = size(avgQL,2);

    %% Number of experiments
    N = experiments;

    %% Number of Gibbs samples
    S = 100;

    mciVariant = 'imci';
    integralRange = [0, max(avgQL(:))];
    thetaStep = integralRange(2) / 400;


    %% assuming uniform prior distribution - TODO add as option
    function p = prior(integralRange, delta)
        p = delta / (integralRange(2) - integralRange(1));
    end

    theta = zeros(S, M, R);
    steps = integralRange(1):thetaStep:integralRange(2);

    for s=1:S
        sampleTheta = theta(s,:,:);
        for i=1:M
            for c=1:R
                logPosteriors = zeros(size(steps));
                for st=1:length(steps)
                    stepTheta = steps(st);
                    logPrior = log(prior(integralRange, thetaStep));
                    gNormalizingConstant = pfqn_mci(squeeze(sampleTheta)', P, Z, N, mciVariant);

                    logPosteriors(st) =  N*avgQL(i, c)*log(stepTheta)-N*log(gNormalizingConstant)+logPrior;
                end

                probs = exp(logPosteriors-max(logPosteriors));
                probs = probs / sum(probs);

                cumulativeProb = cumsum(probs);
                u = rand(1);
                index = find(u<cumulativeProb, 1, 'first');
                sampleTheta(i, c) = steps(index);
            end
        end
        theta(s+1,:,:) = sampleTheta;
    end

    visitPerClass = zeros(R, M);
    for i = 1 : R
      visitPerClass(i, :) = visits{i}(2:M+1, i); %TODO investigate
    end

    theta = arrayfun(@(x) x ./ visitPerClass', theta, 'UniformOutput', false);

    cutoff = round(S / 2);
    theta = theta(cutoff:end);

    thetaAvg = mean(cat(3,theta{:}),3);

    demEst = thetaAvg;
end
