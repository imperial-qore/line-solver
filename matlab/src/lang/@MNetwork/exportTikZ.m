function pdfFile = exportTikZ(self, filePath)
% EXPORTTIKZ Export network diagram to PDF file
%
% PDFFILE = EXPORTTIKZ(FILEPATH) compiles the TikZ code and saves
% the resulting PDF. Returns the File object of the exported PDF.
%
% Requires pdflatex to be installed.
%
% Example:
%   pdfFile = model.exportTikZ('network');
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Convert MATLAB model to Java Network
jnetwork = JLINE.line_to_jline(self);

pdfFile = jnetwork.exportTikZ(filePath);
end
