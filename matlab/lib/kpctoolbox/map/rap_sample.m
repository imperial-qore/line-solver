function sample = rap_sample(RAP, n, xs)
% SAMPLE = RAP_SAMPLE(RAP, N, XS) - Generate random samples from RAP distribution
%
% Generate samples from a Rational Arrival Process (RAP) distribution.
% The marginal distribution of RAP inter-arrival times is a Matrix
% Exponential (ME) distribution, so this function delegates to me_sample.
%
% Input:
%   RAP: RAP distribution as either:
%        - RAP object
%        - Process cell array {H0, H1}
%   N: Number of samples to generate (default: 1)
%   XS: Optional pre-computed grid for CDF evaluation
%       If not provided, auto-generates grid based on mean and variance
%
% Output:
%   SAMPLE: Column vector of N samples from the RAP marginal distribution
%
% Note:
%   This function generates samples from the marginal distribution only.
%   It does not preserve the correlation structure of the RAP.
%   For full RAP simulation with correlations, use map_sample instead.
%
% Examples:
%   rap = RAP([-2, 1; 0.5, -1.5], [0.5, 0.5; 0.5, 0.5]);
%   samples = rap_sample(rap, 10000);
%
%   % Or use process representation directly
%   RAP_proc = {[-2, 1; 0.5, -1.5], [0.5, 0.5; 0.5, 0.5]};
%   samples = rap_sample(RAP_proc, 10000);
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Handle input arguments
if nargin < 2
    n = 1;
end

% Extract process representation
if isa(RAP, 'RAP')
    RAP_proc = RAP.getProcess();
elseif iscell(RAP)
    RAP_proc = RAP;
else
    error('RAP must be either a RAP object or a cell array {H0, H1}');
end

% Delegate to me_sample
% The marginal distribution of RAP is ME with the same {H0, H1} representation
if nargin < 3 || isempty(xs)
    sample = me_sample(RAP_proc, n);
else
    sample = me_sample(RAP_proc, n, xs);
end

end
