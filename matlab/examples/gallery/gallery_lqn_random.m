function model = gallery_lqn_random(seed)
% GALLERY_LQN_RANDOM Randomly generated layered queueing network (reproducible).
% Uses LayeredNetworkGenerator with a fixed default seed so repeated calls
% yield the same model. 1 client, 2 levels, 4 tasks, 2 processors.
if nargin < 1
    seed = 23000;
end
rng(seed);
gen = LayeredNetworkGenerator();
model = gen.generate(1, 2, 4, 2);
end
