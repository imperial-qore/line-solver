function line_error(caller, msg)
%LINE_ERROR Display a plain-text error message with file and line info.
%
%   LINE_ERROR(CALLER, MSG) throws an error with CALLER's name and message,
%   including the source file and line number, in plain text (no hyperlink).

%   Copyright (c) 2012-2026, Imperial College London
%   All rights reserved.

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
end
