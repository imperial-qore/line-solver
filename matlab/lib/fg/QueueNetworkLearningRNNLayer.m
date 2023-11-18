% Adapted from:
% Garbi, G et al. (2020). Learning Queueing Networks by Recurrent Neural Networks
    
classdef QueueNetworkLearningRNNLayer < nnet.layer.Layer % ...
        % & nnet.layer.Formattable ... % (Optional) 
        % & nnet.layer.Acceleratable % (Optional)

    properties
        M
        R
        N
        I
        concurrency
    end

    properties (Learnable)
        % Layer learnable parameters.
        mu
        P
    end

    properties (State)
        hiddenState
    end

    methods
        function layer = QueueNetworkLearningRNNLayer(M,R,concurrency)
            % Create a QueueNetworkLearningRNNLayer

            layer.M = M;
            layer.R = R;
            layer.concurrency = concurrency;
            layer.I = eye(M);

            layer.mu = layer.initializeUniformNonNeg([M,1]);

            P = layer.initializeUniformNonNeg([M,M-1]);
            P = P ./ sum(P, 2);
            oneHot = (find(mod((0:(M^2)-1), M+1)~=0))'==1:M^2;
            layer.P = reshape(reshape(P.', 1, []) * oneHot, M, M);
      
            layer = layer.resetState();
        end

        function parameter = initializeUniformNonNeg(layer, sz)
            a = 0.01;
            b = 10.0;
            parameter = a + (b-a).*rand(sz, 'single');
            parameter = dlarray(parameter);   
        end

        function [Z,state] = predict(layer, X)
            numTimeSteps = size(X,4);
            layer = layer.resetState();
            Z = dlarray(zeros([size(X,1), size(X,2), 1, numTimeSteps]));
            for t=1:numTimeSteps
                currentT = X(1,1,1,t);
                oldT = layer.hiddenState(1);
                deltaT = currentT - oldT;
                pm = abs(layer.mu)'.*(abs(layer.P) - layer.I);
                pred = layer.hiddenState(2:end) + (deltaT*min(layer.hiddenState(2:end), layer.concurrency)) * pm;

                if deltaT == 0
                    pred = X(:,2,1,t)';
                end
                xh_pred = cat(2, [currentT], pred);
                
                layer.hiddenState = xh_pred;
                Z(:,1,:,t) = currentT;
                Z(:,2,:,t) = xh_pred(2:end);
            end
            state = xh_pred;
        end

        function layer = resetState(layer)
            layer.hiddenState = zeros(layer.M+1,1)';
        end
    end
end