function libs = libraries(self, options)
% LIBS = LIBRARIES()
% LIBS = LIBRARIES(OPTIONS)
%
% Third-party libraries this solver will use on its model, as a cell array of
% names, without printing anything. Attribution in LINE is pull-based: nothing
% is written to the console unless you ask, either through this accessor or
% through showLibraryAttribution.
%
% Example:
%   solver = SolverMAM(model);
%   solver.getAvg();
%   libs = solver.libraries()   % {'MAMSolver','BUTools'}
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(options)
    options = self.getOptions;
end
try
    sn = self.getStruct();
catch
    sn = [];
end
libs = self.getLibrariesUsed(sn, options);
end
