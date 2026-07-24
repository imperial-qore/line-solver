pd = pwd;
cwd = fullfile(lineTempDir,'jsim');
cd(cwd)
delete *.jsim
delete *.jsim-result.jsim
cwd = fullfile(lineTempDir,'jmva');
cd(cwd)
delete *.jmva
delete *.jmva-result.jsim
cd(pd)
