function [result, parsed] = getResultsJSIM(self)
% [RESULT, PARSED] = GETRESULTSJSIM()

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

%try
fileName = [getFileName(self),'.jsim-result.jsim'];
filePath = [getFilePath(self),filesep,fileName];
if exist(filePath,'file')
    Pref.Str2Num = 'always';
    parsed = xml_read(filePath,Pref);
else
    line_error(mfilename,'JMT did not output a result file, the simulation has likely failed.');
end
self.result.('solver') = getName(self);
self.result.('model') = parsed.ATTRIBUTE;
self.result.('metric') = {parsed.measure.ATTRIBUTE};

result = self.result;
end
