function line_verbosity(level)
% LINE_VERBOSITY(LEVEL)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

global LINEVerbose;

if nargin<1
    level = VerboseLevel.STD;
end

switch level
    case VerboseLevel.SILENT
        warning OFF
    otherwise
        warning ON BACKTRACE
end
LINEVerbose = level;

end
