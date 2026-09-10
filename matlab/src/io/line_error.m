function line_error(caller, msg, varargin)
%LINE_ERROR Display a plain-text error message with file and line info.
%
%   LINE_ERROR(CALLER, MSG, ...) throws an error with CALLER's name and
%   message, including the source file and line number, in plain text (no
%   hyperlink). Extra arguments are passed to sprintf to format MSG.
%
%   The exception carries NO MATLAB stack, so the console shows the message
%   alone rather than the chain of internal frames that led to it. Every
%   line_error is a diagnostic LINE wrote on purpose and its message already
%   names the throwing function and line, so the frames in between
%   (runAnalyzerChecks -> runAnalyzer -> getAvg -> getAvgTable -> aT ...) are
%   noise to a user who only asked a solver for a result. Genuine MATLAB
%   faults (index out of range, undefined function) are untouched and still
%   report in full.
%
%   Set the verbosity to VerboseLevel.DEBUG to get the stack back, e.g.
%   SolverNC(model,'verbose',VerboseLevel.DEBUG) or
%   GlobalConstants.setVerbose(VerboseLevel.DEBUG). `dbstop if error` stops at
%   the throw site either way, since this is still error().

%   Copyright (c) 2012-2026, Imperial College London
%   All rights reserved.

if coder.target('MATLAB')
    if ~isempty(varargin)
        msg = sprintf(msg, varargin{:});
    end
    msg = strrep(msg, '\n', '');  % Strip out literal '\n' if present
    stack = dbstack;
    if numel(stack) >= 2
        lineNum = stack(2).line;
    else
        lineNum = 1;
    end
    errStr = sprintf('[%s.m @ line %d] %s', caller, lineNum, msg);

    % Defensive: line_error must be able to fire before lineStart has put
    % GlobalConstants on the path or initialized the globals it reads.
    showStack = false;
    try
        verbose = GlobalConstants.getVerbose();
        showStack = ~isempty(verbose) && isscalar(verbose) && verbose >= VerboseLevel.DEBUG;
    catch
    end

    % '%s' rather than errStr as the template: the message is already
    % formatted, and a Windows path in it ("C:\Users\...") would otherwise
    % be read as escape sequences and mangled or rejected.
    if showStack
        error('LINE:Error', '%s', errStr);
    else
        % error() takes the struct's stack verbatim instead of capturing the
        % current one. It must be m-by-1, so a 0-by-0 empty is rejected.
        err.message = errStr;
        err.identifier = 'LINE:Error';
        err.stack = reshape(struct('file',{},'name',{},'line',{}), 0, 1);
        error(err);
    end
else
    error('%s: %s', caller, msg);
end
end
