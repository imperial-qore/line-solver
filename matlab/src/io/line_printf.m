function line_printf(MSG,varargin)
% LINE_PRINTF(MSG, VARARGIN)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% A STATUS ROW IS CLOSED BEFORE ANYTHING ELSE PRINTS. LineStatus rewinds its
% row with backspaces, which walk back over whatever is actually on the
% terminal and cannot cross a newline -- so a warning raised mid-iteration used
% to land inside the row and leave the next rewind short, splicing two
% iterations into one mangled line. Ending the row here gives the warning its
% own line and starts the next iteration on a fresh row. LineStatus writes
% through its own raw path, so this does not recurse.
if LineStatus.isOpen()
    LineStatus.close();
end

if GlobalConstants.Verbose ~= VerboseLevel.SILENT
%     MSG = strrep(MSG,'\n','');
%     MSG = strrep(MSG, '\', '\\');
%     if ~contains(MSG,'...')
%         if contains(MSG,'Summary')
%             fprintf(GlobalConstants.StdOut, sprintf('%s\n',sprintf(MSG, varargin{:})));
%         elseif contains(MSG,'Iter')
%             fprintf(GlobalConstants.StdOut, sprintf('%s',sprintf(MSG, varargin{:})));
%         else
%             fprintf(GlobalConstants.StdOut, '%s\n', sprintf(MSG, varargin{:}));
%         end
%     else
        MSG = sprintf('%s',sprintf(MSG, varargin{:}));
        fprintf(GlobalConstants.StdOut, '%s', MSG);
%    end
end
end