function r = pfqn_cdfun(nvec,cdscaling)
% R = PFQN_CDFUN(NVEC,PHI)

% AMVA-QD class-dependence function

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.
M = size(nvec,1);
r = ones(M,1);
if ~isempty(cdscaling)
    for i = 1:M
        r(i) = 1 / cdscaling{i}(nvec(i,:));
    end
end
end
