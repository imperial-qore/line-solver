function line_verbosity(level)
%LINE_VERBOSITY Sets the global verbosity level for the LINE toolbox.
%
%   LINE_VERBOSITY(LEVEL) sets the verbosity of LINE's output and configures
%   MATLAB's warning behavior accordingly.
%
%   LEVEL should be one of the following (defined in VerboseLevel):
%       - VerboseLevel.SILENT  : Disables warnings and suppresses console output.
%       - VerboseLevel.STD     : Enables standard verbosity and warning backtrace.
%       - VerboseLevel.DEBUG   : Enables debug-level output if supported.
%
%   If LEVEL is not provided, it defaults to VerboseLevel.STD.
%
%   Example:
%       line_verbosity(VerboseLevel.SILENT);
%
%   See also: VerboseLevel, warning

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

global LINEVerbose;

% Default to standard verbosity if no level is provided
if nargin < 1
    level = VerboseLevel.STD;
end

% Adjust MATLAB warning behavior based on selected verbosity level
switch level
    case VerboseLevel.SILENT
        warning off  % Fully suppress warnings
    otherwise
        warning on backtrace  % Enable warnings with backtrace info
end

% Apply selected verbosity level globally
LINEVerbose = level;
end
