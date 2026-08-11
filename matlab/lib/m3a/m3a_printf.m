function m3a_printf(MSG,varargin)
% M3A_PRINTF(MSG, VARARGIN)
%
% Emits an M3A fitting trace message. These messages describe the internal
% progress of the fitting algorithms and are shown only at DEBUG verbosity.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if GlobalConstants.Verbose >= VerboseLevel.DEBUG
    fprintf(GlobalConstants.StdOut, '%s', sprintf(MSG, varargin{:}));
end
end
