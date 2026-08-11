function result = me_open(self, options)
% RESULT = ME_OPEN(OPTIONS)
% Maximum Entropy Method for Open Queueing Networks
%
% Applies the ME algorithm from Kouvatsos (1994) to the model.
% Only supports open queueing networks (no closed classes).
%
% Parameters:
%   options - (optional) struct with fields:
%             .mem_tol     - convergence tolerance (default: 1e-6)
%             .mem_maxiter - maximum iterations (default: 1000)
%             .mem_verbose - print iteration info (default: false)
%
% Returns:
%   result - struct with fields:
%            .QN     - Mean queue lengths [M x R matrix]
%            .UN     - Utilizations [M x R matrix]
%            .RN     - Mean response times [M x R matrix]
%            .TN     - Throughputs (arrival rates) [M x R matrix]
%            .CN     - Completion times [1 x R matrix]
%            .XN     - System throughputs [1 x R matrix]
%
% Reference: Kouvatsos (1994) "Entropy Maximisation and Queueing Network Models"

% Handle optional arguments
if nargin < 2
    options = self.getOptions;
end

% Get model structure
sn = self.getStruct;

% Call solver_nc_mem
[QN, UN, RN, TN, CN, XN, ~] = solver_nc_mem(sn, options);

% Build result structure
result = struct();
result.QN = QN;
result.UN = UN;
result.RN = RN;
result.TN = TN;
result.CN = CN;
result.XN = XN;

end
