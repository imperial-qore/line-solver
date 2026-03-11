function line_error(caller, msg, varargin)
%LINE_ERROR Display a plain-text error message with file and line info.
%
%   LINE_ERROR(CALLER, MSG, ...) throws an error with CALLER's name and
%   message, including the source file and line number, in plain text (no
%   hyperlink). Extra arguments are passed to sprintf to format MSG.

%   Copyright (c) 2012-2026, Imperial College London
%   All rights reserved.

if coder.target('MATLAB')
    if ~isempty(varargin)
        msg = sprintf(msg, varargin{:});
    end
    try
        msg = strrep(msg, '\n', '');  % Strip out literal '\n' if present
        stack = dbstack;

        if numel(stack) >= 2
            lineNum = stack(2).line;
        else
            lineNum = 1;
        end

        filePath = which(caller);
        if isempty(filePath)
            filePath = caller;
        end

        errStr = sprintf('[%s.m @ line %d] %s', caller, lineNum, msg);
        error(errStr);

    catch ME
        throwAsCaller(ME);
    end
else
    error('%s: %s', caller, msg);
end
end
