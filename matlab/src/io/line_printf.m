function line_printf(MSG,varargin)
% LINE_PRINTF(MSG, VARARGIN)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

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