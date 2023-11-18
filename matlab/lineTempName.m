function tmpname = lineTempName(solvername)
if nargin>=1
    tmpname = tempname([lineRootFolder,filesep,'workspace',filesep,solvername]);
else
    tmpname = tempname([lineRootFolder,filesep,'workspace']);
end
if ~exist(tmpname,'dir')
    mkdir(tmpname);
end
end