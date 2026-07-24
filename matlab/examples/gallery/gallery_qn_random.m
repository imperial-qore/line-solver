function model = gallery_qn_random(seed)
% GALLERY_QN_RANDOM Randomly generated mixed queueing network (reproducible).
% Uses NetworkGenerator with a deterministic (cyclic) topology and a fixed
% default seed so repeated calls yield the same model.
% Closed network: 3 queues, 1 delay, 2 closed classes (always stable).
if nargin < 1
    seed = 23000;
end
rng(seed);
gen = NetworkGenerator('schedStrat', 'fcfs', 'routingStrat', 'Probabilities', ...
    'distribution', 'exp', 'cclassJobLoad', 'medium', ...
    'hasVaryingServiceRates', false, 'hasMultiServerQueues', false, ...
    'hasRandomCSNodes', false, 'hasMultiChainCS', false, ...
    'topologyFcn', @NetworkGenerator.cyclicGraph);
model = gen.generate(3, 1, 0, 2);
end
