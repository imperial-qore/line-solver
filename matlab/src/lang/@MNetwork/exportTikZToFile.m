function exportTikZToFile(self, filePath)
% EXPORTTIKZTOFILE Save TikZ/LaTeX source code to a file
%
% EXPORTTIKZTOFILE(FILEPATH) saves the raw TikZ code without compilation,
% useful for manual editing or integration into other LaTeX documents.
%
% Example:
%   model.exportTikZToFile('network.tex')
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Convert MATLAB model to Java Network
jnetwork = JLINE.line_to_jline(self);

jnetwork.exportTikZToFile(filePath);
end
