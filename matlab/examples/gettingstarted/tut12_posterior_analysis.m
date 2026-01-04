% Example 12: Posterior analysis with uncertain parameters
%
% This tutorial demonstrates how to analyze queueing models with parameter
% uncertainty using the Posterior solver. The Posterior solver works with
% Prior distributions that represent uncertainty about model parameters,
% computing posterior distributions of performance metrics.
%
% Scenario: An M/M/1 queue where the service rate is uncertain. We model
% this uncertainty using a Prior distribution with 30 alternatives,
% creating a Gaussian-like distribution of possible service rates.

%% Block 1: Create model with uncertain service rate
model = Network('UncertainServiceModel');

% Create nodes
source = Source(model, 'Source');
queue = Queue(model, 'Queue', SchedStrategy.FCFS);
sink = Sink(model, 'Sink');

% Create job class
jobClass = OpenClass(model, 'Jobs');

% Set arrival rate (lambda = 0.5)
lambda = 0.5;
source.setArrival(jobClass, Exp(lambda));

%% Block 2: Define Prior distribution for uncertain service rate
% Use many alternatives to create a smooth, continuous-looking PDF
% Service rates range from 0.7 to 2.5 with a Gaussian-like prior
numAlternatives = 30;
serviceRates = linspace(0.7, 2.5, numAlternatives);

% Create Gaussian-like prior probabilities centered at mu=1.3
priorMean = 1.3;
priorStd = 0.4;
priorProbs = exp(-0.5 * ((serviceRates - priorMean) / priorStd).^2);
priorProbs = priorProbs / sum(priorProbs);  % Normalize to sum to 1

alternatives = cell(1, numAlternatives);
for i = 1:numAlternatives
    alternatives{i} = Exp(serviceRates(i));
end

prior = Prior(alternatives, priorProbs);
queue.setService(jobClass, prior);

%% Block 3: Complete model topology
model.link(model.serialRouting(source, queue, sink));

fprintf('Model: M/M/1 with uncertain service rate\n');
fprintf('Arrival rate: lambda = %.1f\n', lambda);
fprintf('Number of service rate alternatives: %d\n', numAlternatives);
fprintf('Service rate range: mu in [%.2f, %.2f]\n', min(serviceRates), max(serviceRates));
fprintf('Prior: Gaussian-like centered at mu=%.1f with std=%.1f\n\n', priorMean, priorStd);

%% Block 4: Solve with Posterior wrapper using MVA
post = Posterior(model, @SolverMVA);
post.runAnalyzer();

%% Block 5: Get prior-weighted average results
avgTable = post.getAvgTable()

%% Block 6: Get posterior table with per-alternative results
postTable = post.getPosteriorTable()

%% Block 7: Extract posterior distributions for different metrics
% Get posterior distribution of response time at the queue
respDist = post.getPosteriorDist('R', queue, jobClass);

% Get posterior distribution of queue length at the queue
qlenDist = post.getPosteriorDist('Q', queue, jobClass);

% Get posterior distribution of utilization at the queue
utilDist = post.getPosteriorDist('U', queue, jobClass);

%% Block 8: Extract values and probabilities from the EmpiricalCDF objects
% For response time
respCdfData = respDist.data;  % [CDF, Value]
respValues = respCdfData(:, 2);
respCdf = respCdfData(:, 1);
respProbs = [respCdf(1); diff(respCdf)];  % Convert CDF to PMF

% For queue length
qlenCdfData = qlenDist.data;
qlenValues = qlenCdfData(:, 2);
qlenCdf = qlenCdfData(:, 1);
qlenProbs = [qlenCdf(1); diff(qlenCdf)];

% For utilization
utilCdfData = utilDist.data;
utilValues = utilCdfData(:, 2);
utilCdf = utilCdfData(:, 1);
utilProbs = [utilCdf(1); diff(utilCdf)];

%% Block 9: Print posterior distribution statistics
% Response time statistics
expectedR = sum(respValues .* respProbs);
[~, modeIdxR] = max(respProbs);
fprintf('Response Time (R) at Queue:\n');
fprintf('  Expected Value E[R]: %.4f\n', expectedR);
fprintf('  Mode (most likely): %.4f\n\n', respValues(modeIdxR));

% Queue length statistics
expectedQ = sum(qlenValues .* qlenProbs);
[~, modeIdxQ] = max(qlenProbs);
fprintf('Queue Length (Q) at Queue:\n');
fprintf('  Expected Value E[Q]: %.4f\n', expectedQ);
fprintf('  Mode (most likely): %.4f\n\n', qlenValues(modeIdxQ));

% Utilization statistics
expectedU = sum(utilValues .* utilProbs);
[~, modeIdxU] = max(utilProbs);
fprintf('Utilization (U) at Queue:\n');
fprintf('  Expected Value E[U]: %.4f\n', expectedU);
fprintf('  Mode (most likely): %.4f\n\n', utilValues(modeIdxU));

%% Block 10: Create PDF plots
figure('Name', 'Posterior Distribution PDFs', 'Position', [100, 100, 1200, 400]);

% Plot 1: Response Time PDF
subplot(1, 3, 1);
bar(respValues, respProbs, 1.0, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', [0.1, 0.4, 0.6]);
hold on;
xline(expectedR, 'r-', 'LineWidth', 2.5);
hold off;
xlabel('Response Time (R)');
ylabel('Probability');
title('Posterior PDF: Response Time');
legend({'P(R)', sprintf('E[R] = %.3f', expectedR)}, 'Location', 'best');
grid on;

% Plot 2: Queue Length PDF
subplot(1, 3, 2);
bar(qlenValues, qlenProbs, 1.0, 'FaceColor', [0.4, 0.8, 0.4], 'EdgeColor', [0.2, 0.6, 0.2]);
hold on;
xline(expectedQ, 'r-', 'LineWidth', 2.5);
hold off;
xlabel('Queue Length (Q)');
ylabel('Probability');
title('Posterior PDF: Queue Length');
legend({'P(Q)', sprintf('E[Q] = %.3f', expectedQ)}, 'Location', 'best');
grid on;

% Plot 3: Utilization PDF
subplot(1, 3, 3);
bar(utilValues, utilProbs, 1.0, 'FaceColor', [0.8, 0.4, 0.4], 'EdgeColor', [0.6, 0.2, 0.2]);
hold on;
xline(expectedU, 'r-', 'LineWidth', 2.5);
hold off;
xlabel('Utilization (U)');
ylabel('Probability');
title('Posterior PDF: Utilization');
legend({'P(U)', sprintf('E[U] = %.3f', expectedU)}, 'Location', 'best');
grid on;

sgtitle('Posterior Distribution PDFs for M/M/1 with Uncertain Service Rate');

%% Block 11: Create PDF and CDF plot for response time
figure('Name', 'Response Time: PDF and CDF', 'Position', [150, 150, 1000, 400]);

% Left: Area plot for continuous-looking PDF
subplot(1, 2, 1);
area(respValues, respProbs, 'FaceColor', [0.2, 0.6, 0.8], 'FaceAlpha', 0.7, 'EdgeColor', [0.1, 0.4, 0.6], 'LineWidth', 1.5);
hold on;
xline(expectedR, 'r-', 'LineWidth', 2.5);
[~, modeIdx] = max(respProbs);
plot(respValues(modeIdx), respProbs(modeIdx), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
hold off;
xlabel('Response Time (R)');
ylabel('Probability Density');
title('Posterior PDF');
legend({'f(R)', sprintf('E[R] = %.3f', expectedR), sprintf('Mode = %.3f', respValues(modeIdx))}, 'Location', 'best');
grid on;

% Right: CDF
subplot(1, 2, 2);
plot(respValues, respCdf, 'LineWidth', 2.5, 'Color', [0.2, 0.6, 0.8]);
hold on;
medianIdx = find(respCdf >= 0.5, 1, 'first');
if ~isempty(medianIdx)
    medianR = respValues(medianIdx);
    plot(medianR, 0.5, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    xline(medianR, 'r:', 'LineWidth', 1.5);
end
xline(expectedR, 'g--', 'LineWidth', 2);
yline(0.5, 'k:', 'LineWidth', 1);
hold off;
xlabel('Response Time (R)');
ylabel('Cumulative Probability');
title('Posterior CDF');
if ~isempty(medianIdx)
    legend({'F(R)', sprintf('Median = %.3f', medianR), '', sprintf('E[R] = %.3f', expectedR)}, 'Location', 'best');
else
    legend({'F(R)', sprintf('E[R] = %.3f', expectedR)}, 'Location', 'best');
end
grid on;
xlim([min(respValues)*0.9, max(respValues)*1.1]);
ylim([0, 1.05]);

sgtitle('Response Time Posterior Distribution');
