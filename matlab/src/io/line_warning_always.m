function line_warning_always(caller, MSG, varargin)
% LINE_WARNING_ALWAYS(CALLER, MSG, ...)
%
% Emit a warning without the repeat suppression applied by LINE_WARNING.
%
% LINE_WARNING keeps only the last message and hides an identical repeat for
% 60 seconds. That is right for configuration notices cast once per model, but
% wrong for a warning that reports a correctness limitation of the analysis:
% solving several models in one session would then flag only the first one, and
% the user would read the silence on the others as a clean bill of health.
% Warnings that say "these numbers are not exact" must be raised for every model
% they apply to, so they go through here instead. Verbosity gating is unchanged:
% SILENT still silences everything.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~coder.target('MATLAB')
    return;  % No-op in codegen mode
end

if GlobalConstants.Verbose == VerboseLevel.SILENT
    return
end

errmsg = sprintf(MSG, varargin{:});
% '%s\n' rather than passing the assembled string as the format itself: the
% message carries user-supplied names and must not be reinterpreted as a
% format specification.
line_printf('%s\n', sprintf('Warning [%s.m]: %s', caller, errmsg));
end
