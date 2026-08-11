function line_debug(varargin)
% LINE_DEBUG(MSG, VARARGIN)
% LINE_DEBUG(OPTIONS, MSG, VARARGIN)
%
% Print debug message if verbose level is DEBUG.
% If OPTIONS struct is passed as first argument, checks options.verbose.
% Otherwise checks GlobalConstants.Verbose.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Check if first argument is an options struct with verbose field
if nargin >= 1 && isstruct(varargin{1}) && isfield(varargin{1}, 'verbose')
    options = varargin{1};
    MSG = varargin{2};
    args = varargin(3:end);
    isDebug = (options.verbose == VerboseLevel.DEBUG) || ...
              (GlobalConstants.Verbose == VerboseLevel.DEBUG);
else
    MSG = varargin{1};
    args = varargin(2:end);
    isDebug = (GlobalConstants.Verbose == VerboseLevel.DEBUG);
end

if isDebug
    MSG = sprintf('[DEBUG] %s', sprintf(MSG, args{:}));
    line_printf('%s\n', MSG);
end
end