function result = qsys_mapd1(D0, D1, s, varargin)
% QSYS_MAPD1 Analyzes a MAP/D/1 queue using Q-MAM.
%
% RESULT = QSYS_MAPD1(D0, D1, S) analyzes a MAP/D/1 queue with:
%   D0 - MAP hidden transition matrix (n x n)
%   D1 - MAP arrival transition matrix (n x n)
%   S  - Deterministic service time (positive scalar)
%
% RESULT = QSYS_MAPD1(..., 'maxNumComp', N) sets max queue length components (default 1000)
% RESULT = QSYS_MAPD1(..., 'numSteps', K) sets waiting time distribution granularity (default 1)
%
% This is a convenience wrapper for qsys_mapdc with c=1.
%
% See also qsys_mapdc, Q_CT_MAP_D_C, qsys_mapm1

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

result = qsys_mapdc(D0, D1, s, 1, varargin{:});

end
