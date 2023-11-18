function tmpdir = lineTempDir
    tmpdir = [lineRootFolder,filesep,'workspace',filesep];
    if ~exist(tmpdir,'dir')
        mkdir(tmpdir);
    end
end