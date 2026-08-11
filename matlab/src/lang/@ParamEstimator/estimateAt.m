function estVal = estimateAt(self, nodes)

sn = self.model.getStruct;

if ~iscell(nodes)
    nodes = {nodes};
end

switch self.options.method
    case 'ubr' % utilization based regression
        estVal = estimator_ubr(self, nodes);
    case 'ubo' % utilization based optimization
        estVal = estimator_ubo(self, nodes);
    case 'erps' % extended regression for PS
        estVal = estimator_erps(self, nodes);
    case 'ekf' % Extended Kalman Filter Estimation
        estVal = estimator_ekf(self, nodes);
    case 'mcmc' % Gibbs Sampling MCMC
        estVal = estimator_mcmc(self, nodes);
    case 'mle' % Maximum Likelihood Estimation
        estVal = estimator_mle(self, nodes);
    case 'rnn' % Explainable RNN Estimation
        estVal = estimator_rnn(self, nodes);
    case 'mlps' % Maximum Likelihood for PS
        estVal = estimator_mlps(self, nodes);
    case 'fmlps' % Fluid Maximum Likelihood for PS
        estVal = estimator_fmlps(self, nodes);
    case 'qmle' % Quick MLE
        estVal = estimator_qmle(self, nodes);
    case 'gibbs' % Gibbs Sampling
        estVal = estimator_gibbs(self, nodes);
    otherwise
        error('Unknown inference method: %s.', self.options.method);
end

% update the model parameters
for n=1:size(nodes, 2)
    svcProc = nodes{n}.getService;
    for r=1:sn.nclasses
        svcProc{r}.setMean(estVal(n, r));
    end
end
self.model.reset;

end
