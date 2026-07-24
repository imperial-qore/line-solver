function out = getJMVATempPath(self)
% OUT = GETJMVATEMPPATH()

if isempty(self.filePath) || isempty(self.fileName)
    self.filePath = lineTempName('jmva');
    self.fileName = 'model';
end
fname = [self.fileName,'.jmva'];
out = [self.filePath,filesep,fname];
end
