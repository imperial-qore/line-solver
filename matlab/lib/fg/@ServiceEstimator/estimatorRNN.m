function estVal = estimatorRNN(self, nodes)
    % Explainable RNN techinque to learn queue structure and parameters
    % with the aim of service demand estimation
    %
    % Adapted using:
    % Garbi, G et al. (2020). Learning Queueing Networks by Recurrent Neural Networks
    % Copyright (c) 2012-2024, Imperial College London
    % All rights reserved.
    % This code is released under the 3-Clause BSD License.
    
    sn = self.model.getStruct;

    % use whole model nodes for RNN estimation
    nodes = self.model.getNodes();
    qLenTs = {};
    qLenTrace = {};
    % obtain per class metrics
    minSampleCount = inf;
    for n=1:size(nodes,1)
        node = nodes{n};
        numServers(n) = node.getNumberOfServers();
        for r=1:sn.nclasses
            samples = self.getQLen(node,self.model.classes{r}); % queue-length
            if isempty(samples)
                error('Queue-length data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
            else
                for d=1:length(samples)
                    qLen = samples{d};

                    if length(qLenTs) < d
                        qLenTs{d} = {};
                        qLenTrace{d} = {};
                    end

                    qLenTs{d}{n,r} = qLen.t;
                    qLenTrace{d}{n,r} = qLen.data;

                    if size(qLen.data, 1) < minSampleCount
                        minSampleCount = size(qLen.data, 1);
                    end
                end
            end
        end
    end

    traces = [];
    for n=1:length(qLenTs)
        traceQL = cell2mat(qLenTrace{n});
        traceQL = reshape(traceQL, size(traceQL,1)/size(qLenTrace{n}, 1), size(qLenTrace{n}, 1), sn.nclasses);
        traceTs = cell2mat(qLenTs{n});
        traceTs = reshape(traceTs, size(traceTs,1)/size(qLenTs{n}, 1), size(qLenTs{n}, 1), sn.nclasses);
        traceQL = cat(3, traceTs, traceQL);
        traces(end+1, :, :, :) = traceQL(1:minSampleCount,:,:);
    end

    estVal = rnn_data(traces, numServers);
end

% rnn procedure
function demEst = rnn_data(avgQL, numServers)

    %% number of resources
    M = size(avgQL, 3);

    %% number of classes
    R = size(avgQL, 4)-1;
    
    %% Number of samples
    S = size(avgQL, 2);

    %% Number of experiments
    traceCount = size(avgQL, 1);
    

    numEpochs = 2;
    numIterationsPerEpoch = 50;
   
    qrnnLayer = QueueNetworkLearningRNNLayer(M,R,numServers);
    layers = [
        sequenceInputLayer([M, R+1, 1]);
        qrnnLayer];

    net = dlnetwork(layers, initialize=true);

    % Define max absolute percentage error loss function
    function [loss, gradients, state] = modelLoss(net, X, T)
        % Forward data through the dlnetwork object.
        [Y, state] = forward(net,X);
        
        predErr = abs(T(:, 2:end, 1, :) - Y(:, 2:end, 1, :));

        N = mean(sum(X(:, 2:end, 1, :), 1));
        maxErr = max(sum(predErr,1) ./ (2.*N), [], 4);
        loss = 100*maxErr;
    
        % Compute gradients.
        gradients = dlgradient(loss, net.Learnables);
    end
    lr = 0.1;
    
    % Loop over epochs.
    for epoch = 1:numEpochs
         iteration = 1;
         trailingAvg  = [];
         trailingAverageSq = [];
        % Loop over mini-batches.
        for i = 1:numIterationsPerEpoch
            iteration = iteration + 1;
            
            qrnnLayer.resetState();

            exp = randi(traceCount);
            trace = squeeze(avgQL(exp, :, :, :));

            X = dlarray(trace,'TSSC');
            T = dlarray(trace,'TSSC');
            % Evaluate model loss and gradients.
            [loss,gradients,state] = dlfeval(@modelLoss,net,X,T);
            net.State = state;
            disp(loss);

            % Update learnable parameters.
            [net,trailingAvg,trailingAverageSq] = adamupdate(net, gradients, trailingAvg,trailingAverageSq,iteration, lr);
    
        end
    end

    mu = net.Learnables.Value(1);
    demEst = extractdata(mu{1});
end
