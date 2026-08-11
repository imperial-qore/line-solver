function estVal = estimator_rnn(self, nodes)
    % Explainable RNN technique to learn queue structure and parameters
    % with the aim of service demand estimation.
    %
    % Delegates training to a PyTorch-based RNN via a Python bridge script,
    % removing the dependency on the MATLAB Deep Learning Toolbox.
    %
    % Adapted from:
    % Garbi, G et al. (2020). Learning Queueing Networks by Recurrent Neural Networks
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    % This code is released under the 3-Clause BSD License.

    sn = self.model.getStruct;

    % Use whole model nodes for RNN estimation
    nodes = self.model.getNodes();
    qLenTs = {};
    qLenTrace = {};
    % Obtain per-class metrics
    minSampleCount = inf;
    for n = 1:size(nodes, 1)
        node = nodes{n};
        numServers(n) = node.getNumberOfServers();
        for r = 1:sn.nclasses
            samples = self.getQLen(node, self.model.classes{r});
            if isempty(samples)
                error('Queue-length data for node %d in class %d is missing.', self.model.getNodeIndex(node), r);
            else
                for d = 1:length(samples)
                    qLen = samples{d};

                    if length(qLenTs) < d
                        qLenTs{d} = {};
                        qLenTrace{d} = {};
                    end

                    qLenTs{d}{n, r} = qLen.t;
                    qLenTrace{d}{n, r} = qLen.data;

                    if size(qLen.data, 1) < minSampleCount
                        minSampleCount = size(qLen.data, 1);
                    end
                end
            end
        end
    end

    % Build 4D traces array: (traceCount x S x M x R+1)
    traces = [];
    for n = 1:length(qLenTs)
        traceQL = cell2mat(qLenTrace{n});
        traceQL = reshape(traceQL, size(traceQL, 1) / size(qLenTrace{n}, 1), size(qLenTrace{n}, 1), sn.nclasses);
        traceTs = cell2mat(qLenTs{n});
        traceTs = reshape(traceTs, size(traceTs, 1) / size(qLenTs{n}, 1), size(qLenTs{n}, 1), sn.nclasses);
        traceQL = cat(3, traceTs, traceQL);
        traces(end + 1, :, :, :) = traceQL(1:minSampleCount, :, :);
    end

    estVal = rnn_data(traces, numServers);
end

% Call Python PyTorch RNN via bridge script
function demEst = rnn_data(avgQL, numServers)
    % Locate the Python bridge script relative to this file
    thisDir = fileparts(mfilename('fullpath'));
    repoRoot = fullfile(thisDir, '..', '..', '..', '..');
    bridgeScript = fullfile(repoRoot, 'python', 'line_inference', 'api', 'rnn_bridge.py');

    % Create temporary files for data exchange
    inputFile = [tempname, '.mat'];
    outputFile = [tempname, '.mat'];

    % Save input data
    traces = avgQL; %#ok<NASGU>
    save(inputFile, 'traces', 'numServers', '-v7');

    % Call Python bridge
    pythonCmd = sprintf('python3 %s %s %s', bridgeScript, inputFile, outputFile);

    % Add the Python package to PYTHONPATH
    pythonPath = fullfile(repoRoot, 'python');
    envCmd = sprintf('PYTHONPATH=%s:$PYTHONPATH %s', pythonPath, pythonCmd);

    [status, cmdout] = system(envCmd);

    % Clean up input file
    delete(inputFile);

    if status ~= 0
        if exist(outputFile, 'file')
            delete(outputFile);
        end
        error('RNN Python bridge failed (exit code %d):\n%s', status, cmdout);
    end

    % Load results
    result = load(outputFile);
    delete(outputFile);

    demEst = result.demandEst(:)';
end
