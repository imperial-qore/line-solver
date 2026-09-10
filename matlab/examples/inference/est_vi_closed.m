% Variational inference for Markovian queueing networks, on a closed loop.
%
% Method 'vi' (I. Perez, G. Casale, Adv. Appl. Prob. 53(3), 2021) infers
% service rates from NOISY QUEUE-LENGTH READINGS taken over time: each reading
% is exact with probability 1-epsilon and uniform over the remaining feasible
% values otherwise. Unlike the other estimators it returns a conjugate Gamma
% POSTERIOR per rate, not only a point estimate.
%
% It reads the queue lengths of EVERY station, not only the estimated one: the
% transition counts the method is written in are pinned by the whole picture.
clear node jobclass solver AvgTable

%% define model, with the queue rate to be estimated
N = 10;
model = Network('model');
node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
jobclass{1} = ClosedClass(model, 'Class1', N, node{1}, 0);
node{1}.setService(jobclass{1}, Exp(0.5));
node{2}.setService(jobclass{1}, Exp(2.0));   % starting point of the estimate
model.link(Network.serialRouting(node{1}, node{2}));

%% queue-length readings, one per unit time, 10%% of them faulty
ts = (1:10)';
qlen = [1 2 3 2 4 3 5 4 3 4]';

%% estimate the service rate of the queue
options = ParamEstimator.defaultOptions;
options.method = 'vi';
options.epsilon = 0.1;      % probability that a reading is faulty
options.prior_shape = 2;    % Gamma prior shape; the rate is set from the model
options.ngrid = 51;         % time grid of the backward and forward passes
options.nsamples = 32;      % lattice points per marginal
options.ymax = 60;          % transition-count truncation
options.iter_max = 5;
se = ParamEstimator(model, options);
se.addSamples(SampledMetric(MetricType.QLen, ts, N - qlen, node{1}, jobclass{1}));
se.addSamples(SampledMetric(MetricType.QLen, ts, qlen, node{2}, jobclass{1}));
estVal = se.estimateAt(node{2})

fprintf(1, 'posterior service rate ~ Gamma(%.4f, %.4f), mean %.4f\n', ...
    se.options.posterior(1,1), se.options.posterior(1,2), ...
    se.options.posterior(1,1)/se.options.posterior(1,2));
fprintf(1, 'evidence lower bound over the iterations: %s\n', ...
    mat2str(round(se.options.bound', 3)));

%% solve the model the estimate has been written into
solver = {};
solver{end+1} = MVA(model);
AvgTable = cell(1,length(solver));
for s=1:length(solver)
    fprintf(1,'SOLVER: %s\n',solver{s}.getName());
    AvgTable{s} = solver{s}.getAvgTable();
    AvgTable{s}
end
