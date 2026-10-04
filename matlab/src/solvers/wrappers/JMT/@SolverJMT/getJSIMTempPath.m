function out = getJSIMTempPath(self)
% OUT = GETJSIMTEMPPATH()

if isempty(self.filePath) || isempty(self.fileName)
    self.filePath = lineTempName('jsim');
    self.fileName = 'model';
end
fname = [self.fileName,'.jsim'];
out = [self.filePath,filesep,fname];
end