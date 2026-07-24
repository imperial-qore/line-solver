function vals = matlab_ilt_matrix(fun, T, maxFnEvals, method)
%MATLAB_ILT_MATRIX  Numerical inverse Laplace transform of a matrix-valued function.
%   VALS = MATLAB_ILT_MATRIX(FUN, T, MAXFNEVALS, METHOD) inverts the
%   matrix-valued Laplace transform FUN (a handle s -> matrix) at the time
%   points in T, using the Abate-Whitt framework with the same eta/beta weights
%   as the scalar MATLAB_ILT. This variant evaluates FUN once per node and
%   accumulates the full matrix, avoiding one scalar inversion per entry.
%
%   METHOD is one of 'cme' (default), 'euler', or 'gaver'. Returns a
%   numel(T) x nr x nc array where [nr,nc] = size(FUN(.)).
%
%   The CME parameter table is loaded from iltcme.json in this folder (shared
%   with MATLAB_ILT).
if nargin < 4 || isempty(method); method = 'cme'; end

persistent cmeParamsMat cmeParamsMatPath
here = fileparts(mfilename('fullpath'));
jsonPath = fullfile(here, 'iltcme.json');

switch lower(method)
    case 'cme'
        if isempty(cmeParamsMat) || ~strcmp(cmeParamsMatPath, jsonPath)
            cmeParamsMat = jsondecode(fileread(jsonPath));
            cmeParamsMatPath = jsonPath;
        end
        params = cmeParamsMat(1);
        for i = 2:length(cmeParamsMat)
            if cmeParamsMat(i).cv2 < params.cv2 && (cmeParamsMat(i).n + 1) <= maxFnEvals
                params = cmeParamsMat(i);
            end
        end
        a = params.a(:);
        b = params.b(:);
        c = params.c;
        mu1 = params.mu1;
        omega = params.omega;
        nn = params.n;
        eta  = [c * mu1; (a + 1i*b) * mu1];
        beta = [1;        1 + 1i * omega * (1:nn).'] * mu1;
    case 'euler'
        n_euler = floor((maxFnEvals-1)/2);
        eta = [0.5, ones(1, n_euler), zeros(1, n_euler-1), 2^-n_euler];
        for k = 1:n_euler-1
            eta(2*n_euler-k + 1) = eta(2*n_euler-k + 2) + ...
                exp(sum(log(1:n_euler)) - n_euler*log(2) - sum(log(1:k)) - sum(log(1:(n_euler-k))));
        end
        kidx = 0:2*n_euler;
        beta = n_euler*log(10)/3 + 1i*pi*kidx;
        eta  = (10^((n_euler)/3))*(1-mod(kidx, 2)*2) .* eta;
        eta = eta(:); beta = beta(:);
    case 'gaver'
        if mod(maxFnEvals,2)==1
            maxFnEvals = maxFnEvals - 1;
        end
        ndiv2 = maxFnEvals/2;
        eta = zeros(maxFnEvals,1);
        beta = zeros(maxFnEvals,1);
        for k = 1:maxFnEvals
            inside_sum = 0.0;
            for j = floor((k+1)/2):min(k,ndiv2)
                inside_sum = inside_sum + exp((ndiv2+1)*log(j) - sum(log(1:(ndiv2-j))) + ...
                    sum(log(1:2*j)) - 2*sum(log(1:j)) - sum(log(1:(k-j))) - sum(log(1:(2*j-k))));
            end
            eta(k) = log(2.0)*(-1)^(k+ndiv2)*inside_sum;
            beta(k) = k * log(2.0);
        end
    otherwise
        line_error(mfilename, sprintf('Unknown inverse Laplace method "%s". Supported: cme, euler, gaver', method));
end

% Probe output size with one cheap evaluation.
x0 = T(1);
M0 = fun(beta(1) / x0);
[nr, nc] = size(M0);

vals = zeros(numel(T), nr, nc);
nBeta = numel(eta);
for ii = 1:numel(T)
    x = T(ii);
    acc = zeros(nr, nc);
    for k = 1:nBeta
        acc = acc + eta(k) * fun(beta(k) / x);
    end
    vals(ii, :, :) = real(acc) / x;
end
end
