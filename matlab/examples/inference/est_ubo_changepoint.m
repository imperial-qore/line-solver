clear node jobclass
%% Example: UBO sliding-window change-point detection (open network)
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
jobclass{2} = OpenClass(model, 'Class2', 0);

node{1}.setService(jobclass{1}, Exp.fitMean(1.0));
node{1}.setService(jobclass{2}, Exp.fitMean(1.0));
node{2}.setService(jobclass{1}, Exp(NaN));
node{2}.setService(jobclass{2}, Exp(NaN));
node{3}.setArrival(jobclass{1}, Exp(1.0));
node{3}.setArrival(jobclass{2}, Exp(0.5));

P = model.initRoutingMatrix;
P{1,1} = [0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0];
P{2,2} = [0,1,0,0; 0,0,0,1; 1,0,0,0; 0,0,0,0];
model.link(P);

%% Generate data with change point
% Class 1 demand changes: D1=0.1 -> 0.3
% Class 2 demand stays: D2=0.2
n = 200;
changeT = 100;
D1_before = 0.1; D1_after = 0.3;
D2 = 0.2;
lambda1 = 1.0; lambda2 = 0.5;
noise = 0.01;

ts = (1:n)';
arvr1_samples = lambda1*ones(n,1) + randn(n,1)*noise;
arvr2_samples = lambda2*ones(n,1) + randn(n,1)*noise;

util_samples = zeros(n,1);
respt1_samples = zeros(n,1);
respt2_samples = zeros(n,1);
for t = 1:n
    if t <= changeT
        d1 = D1_before;
    else
        d1 = D1_after;
    end
    U = lambda1*d1 + lambda2*D2 + randn(1)*noise;
    U = max(0.01, min(0.95, U));
    util_samples(t) = U;
    respt1_samples(t) = d1 / (1 - U) + randn(1)*noise*0.1;
    respt2_samples(t) = D2 / (1 - U) + randn(1)*noise*0.1;
end

%% Sliding-window estimation
W = 30;
est1 = zeros(n - W + 1, 1);
est2 = zeros(n - W + 1, 1);

for t = W:n
    idx = (t-W+1):t;

    node{2}.setService(jobclass{1}, Exp(NaN));
    node{2}.setService(jobclass{2}, Exp(NaN));
    model.reset();

    estoptions = ParamEstimator.defaultOptions;
    estoptions.method = 'ubo';
    se = ParamEstimator(model, estoptions);

    se.addSamples(SampledMetric(MetricType.ArvR, ts(idx), arvr1_samples(idx), node{2}, jobclass{1}));
    se.addSamples(SampledMetric(MetricType.ArvR, ts(idx), arvr2_samples(idx), node{2}, jobclass{2}));
    se.addSamples(SampledMetric(MetricType.RespT, ts(idx), respt1_samples(idx), node{2}, jobclass{1}));
    se.addSamples(SampledMetric(MetricType.RespT, ts(idx), respt2_samples(idx), node{2}, jobclass{2}));
    se.addSamples(SampledMetric(MetricType.Util, ts(idx), util_samples(idx), node{2}));
    se.interpolate();
    estVal = se.estimateAt(node{2});
    est1(t - W + 1) = estVal(1);
    est2(t - W + 1) = estVal(2);
end

%% Report results
fprintf(1, '\n=== UBO Change-Point Detection (2-class) ===\n');
fprintf(1, 'Class 1: D=%.1f (t<=100), D=%.1f (t>100)\n', D1_before, D1_after);
fprintf(1, 'Class 2: D=%.1f (constant)\n', D2);
fprintf(1, 'Window size: %d\n\n', W);

checkpoints = [W, 50, 80, 100, 110, 120, 130, 150, 180, 200];
fprintf(1, '%6s  %8s  %8s  %8s  %8s\n', 't', 'Est D1', 'True D1', 'Est D2', 'True D2');
fprintf(1, '%6s  %8s  %8s  %8s  %8s\n', '------', '--------', '--------', '--------', '--------');
for i = 1:length(checkpoints)
    t = checkpoints(i);
    if t >= W && t <= n
        e1 = est1(t - W + 1);
        e2 = est2(t - W + 1);
        if t <= changeT
            trueD1 = D1_before;
        else
            trueD1 = D1_after;
        end
        fprintf(1, '%6d  %8.4f  %8.4f  %8.4f  %8.4f\n', t, e1, trueD1, e2, D2);
    end
end

%% Plot
figure;
tAxis = (W:n)';
subplot(2,1,1);
plot(tAxis, est1, 'b-', 'LineWidth', 1.5); hold on;
yline(D1_before, 'r--', 'LineWidth', 1);
yline(D1_after, 'r--', 'LineWidth', 1);
xline(changeT, 'k:', 'Change point', 'LineWidth', 1);
ylabel('Demand'); title('Class 1 (changes D=0.1 -> D=0.3)');
grid on;

subplot(2,1,2);
plot(tAxis, est2, 'b-', 'LineWidth', 1.5); hold on;
yline(D2, 'r--', 'LineWidth', 1);
xline(changeT, 'k:', 'Change point', 'LineWidth', 1);
xlabel('Time'); ylabel('Demand'); title('Class 2 (constant D=0.2)');
grid on;
