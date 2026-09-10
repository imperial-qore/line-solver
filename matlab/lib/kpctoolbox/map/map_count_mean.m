function [m] = map_count_mean(map,t)
% Computes the mean of the counting process, at resolution t, for the
% given MAP.
% Input:
% - mmap: Markovian Arrival Process
% - t: the period considered for each sample of the counting process
% Output:
% - m: mean arrivals in (0,t]

% arrival rate
l = map_lambda(map);

% mean
m = l * t(:);
end