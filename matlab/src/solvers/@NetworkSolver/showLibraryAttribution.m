function showLibraryAttribution(self, options)
% SHOWLIBRARYATTRIBUTION(OPTIONS)
%
% Print the third-party libraries the solver will use, at most once per
% session. This is not called automatically: attribution is pull-based, so a
% user asks for it (see also getLibrariesUsed, which returns the list without
% printing).
% getLibrariesUsed is a static method of the concrete solver class and is
% resolved through the instance, so the block no longer has to be copied into
% each runAnalyzer just to name that class.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2 || isempty(options)
    options = self.getOptions;
end
if options.verbose == VerboseLevel.SILENT || GlobalConstants.isLibraryAttributionShown()
    return
end
try
    sn = self.getStruct();
catch
    sn = [];
end
libs = self.getLibrariesUsed(sn, options);
if ~isempty(libs)
    line_printf('The solver will leverage %s.\n', strjoin(libs, ', '));
    GlobalConstants.setLibraryAttributionShown(true);
end
end
