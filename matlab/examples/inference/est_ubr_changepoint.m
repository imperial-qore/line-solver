clear node jobclass
%% Example: UBR sliding-window change-point detection (open network)
%
% Service demand changes from D=0.3 to D=0.6 at t=100.
% A sliding window of size W is used to estimate demands online.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%% define model
model = Network('model');

node{1} = Delay(model, 'Delay');
node{2} = Queue(model, 'Queue1', SchedStrategy.FCFS);
node{3} = Source(model, 'Source');
node{4} = Sink(model, 'Sink');

jobclass{1} = OpenClass(model, 'Class1', 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{2}.setService(jobclass{1}, Exp(NaN));
node{3}.setArrival(jobclass{1}, Exp(1.0));

P = model.initRoutingMatrix;
P{1,1} = [0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0];
model.link(P);

%% Generate data with change point
n = 200;
changeT = 100;
D1 = 0.3;  % demand before change
D2 = 0.6;  % demand after change
lambda = 1.0;
noise = 0.02;

ts = (1:n)';
arvr_samples = lambda*ones(n,1) + randn(n,1)*noise;

% Utilization and response time change with demand
util_samples = zeros(n,1);
respt_samples = zeros(n,1);
for t = 1:n
    if t <= changeT
        D = D1;
    else
        D = D2;
    end
    U = lambda * D + randn(1)*noise;
    U = max(0.01, min(0.95, U));
    util_samples(t) = U;
    respt_samples(t) = D / (1 - U) + randn(1)*noise*0.1;
end

%% Sliding-window estimation
W = 30; % window size
estimates = zeros(n - W + 1, 1);

for t = W:n
    idx = (t-W+1):t;

    node{2}.setService(jobclass{1}, Exp(NaN));
    model.reset();

    estoptions = ParamEstimator.defaultOptions;
    estoptions.method = 'ubr';
    se = ParamEstimator(model, estoptions);

    se.addSamples(SampledMetric(MetricType.ArvR, ts(idx), arvr_samples(idx), node{2}, jobclass{1}));
    se.addSamples(SampledMetric(MetricType.Util, ts(idx), util_samples(idx), node{2}));
    se.interpolate();
    estVal = se.estimateAt(node{2});
    estimates(t - W + 1) = estVal(1);
end

%% Report results
tAxis = (W:n)';
fprintf(1, '\n=== UBR Change-Point Detection ===\n');
fprintf(1, 'True demand: D=%.1f (t<=100), D=%.1f (t>100)\n', D1, D2);
fprintf(1, 'Window size: %d\n\n', W);

% Report estimates at key time points
checkpoints = [W, 50, 80, 100, 110, 120, 130, 150, 180, 200];
fprintf(1, '%6s  %10s  %10s\n', 't', 'Estimated', 'True');
fprintf(1, '%6s  %10s  %10s\n', '------', '----------', '----------');
for i = 1:length(checkpoints)
    t = checkpoints(i);
    if t >= W && t <= n
        est = estimates(t - W + 1);
        if t <= changeT
            trueD = D1;
        else
            trueD = D2;
        end
        fprintf(1, '%6d  %10.4f  %10.4f\n', t, est, trueD);
    end
end

%% Plot
figure;
plot(tAxis, estimates, 'b-', 'LineWidth', 1.5); hold on;
yline(D1, 'r--', sprintf('D=%.1f', D1), 'LineWidth', 1);
yline(D2, 'r--', sprintf('D=%.1f', D2), 'LineWidth', 1);
xline(changeT, 'k:', 'Change point', 'LineWidth', 1);
xlabel('Time');
ylabel('Estimated Demand');
title('UBR Sliding-Window Change-Point Detection');
legend('UBR estimate', 'Location', 'southeast');
grid on;
