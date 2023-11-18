function estVal = estimateAt(self, nodes)

sn = self.model.getStruct;

switch self.options.method
    case 'ubr' % utilization based regression
        estVal = estimatorUBR(self, nodes);
    case 'ubo' % utilization based optimization
        estVal = estimatorUBO(self, nodes);
    case 'ubo2' % altternative utilization based optimization
        estVal = estimatorUBO2(self, nodes);
    case 'ubo3' % utilization based optimization - multi resource attempt
        estVal = estimatorUBO3(self, nodes);
    case 'erps' % extended regression for PS
        estVal = estimatorERPS(self, nodes);
    case 'ekf' % Extended Kalman Filter Estimation
        estVal = estimatorEKF(self, nodes);
    case 'mcmc' % Gibbs Sampling MCMC
        estVal = estimatorMCMC(self, nodes);
    case 'mle' % Maximum Likelihood Estimation
        estVal = estimatorMLE(self, nodes);
    case 'rnn' % Explainable RNN Estimation
        estVal = estimatorRNN(self, nodes);
    otherwise
        error('Unknown inference method: %s.', self.options.method);
end

if ~iscell(nodes)
    nodes = {nodes};

% update the model parameters
for n=1:size(nodes, 2)
    svcProc = nodes{n}.getService;
    for r=1:sn.nclasses
        svcProc{r}.updateMean(estVal(n, r));
    end
end

end
