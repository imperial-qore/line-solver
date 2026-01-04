function fileName = getFileName(self)
% GETFILENAME Get the output file name
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if isempty(self.fileName)
    fileName = 'model';
else
    fileName = self.fileName;
end
end
