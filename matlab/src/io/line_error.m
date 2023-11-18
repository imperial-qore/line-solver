function line_error(caller,MSG)
% LINE_ERROR(CALLER, ERRMSG)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.
try
    MSG = strrep(MSG,'\n','');
    x=dbstack;
    error('[<a href="matlab:opentoline(''%s'',%d)">%s.m@%d</a>] %s',which(caller),x(2).line,caller,x(2).line,MSG);
catch ME
    throwAsCaller(ME)
end
end
