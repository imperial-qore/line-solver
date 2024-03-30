function [sts]=dtmc_simulate(P, pi0, numSamples)
% [sts]=dtmc_simulate(P, pi0, n)

% Ensure PMF sums to 1
pi0 = pi0 / sum(pi0);
initialState = find(cumsum(pi0) >= rand, 1, 'first');
sts = zeros(1, numSamples); % Preallocate the array of states
sts(1) = initialState; % Set the initial state

CP = cumsum(P,2);
for i = 2:numSamples
    sts(i) = find(CP(sts(i-1),:) >= rand, 1); % Update the state sequence
end
end