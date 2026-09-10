function RD = getTranCdfRespT(self, R)
% RD = GETTRANCDFRESPT(R)
%
% Empirical response time CDF under the transient getter's name.
%
% The LDES engine keeps ONE set of per-job response time samples, so this is
% the same measured ecdf getCdfRespT reports -- the C++ CLI serves the same
% curve under -a cdf, -a tran-cdf-respt and -a tran-cdf-passt, and the JAR's
% getTranCdfRespT delegates the same way. RD is an (nstations x nclasses)
% cell of [F(t), t] matrices.
%
% A LayeredNetwork is refused: the layered engine reports steady state only
% and its response times belong to entries, not to a (station, class) grid.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isa(self.model, 'LayeredNetwork')
    line_error(mfilename, ['getTranCdfRespT is indexed by (station, class) and this solver holds a ' ...
        'LayeredNetwork, whose response times belong to entries and whose engine reports ' ...
        'steady state only. Use getCdfRespTLN.']);
end

if nargin < 2
    RD = self.getCdfRespT();
else
    RD = self.getCdfRespT(R);
end
end
