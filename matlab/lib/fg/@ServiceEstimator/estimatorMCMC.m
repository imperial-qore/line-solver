function estVal = estimatorMCMC(self, nodes)
    % Gibbs Sampling MCMC-based optimization
    %
    % Copyright (c) 2012-2022, Imperial College London
    % All rights reserved.
    % This code is released under the 3-Clause BSD License.
    
    sn = self.model.getStruct;
    
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
    
    estVal = mcmc_data(avgQL, self.model.getStruct.visits, experiments, self.options.iter_max);
end

% mcmc procedure based on the comon data format
function demEst = mcmc_data(avgQL, visits, experiments, ITERMAX)

    %% number of resources
    M = size(avgQL,1);
    %% number of classes
    R = size(avgQL,2);
    %% population counts
    P = ones(1,R);
    %% think times
    Z = ones(1,R);
    
    %% Number of experiments
    N = experiments;
    
    %% Number of Gibbs samples
    S = 100;

    mciVariant = 'imci';
    integralRange = [0, 0.2];
    thetaStep = 0.0005;
    
    
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
