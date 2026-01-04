workspaceFolder = [lineRootFolder, filesep, 'workspace/'];
files = dir(workspaceFolder);
for k = 1:length(files)
    if ~files(k).isdir || ~ismember(files(k).name, {'.', '..'})
        if isfolder(fullfile(workspaceFolder, files(k).name))
            rmdir(fullfile(workspaceFolder, files(k).name), 's');
        else
            delete(fullfile(workspaceFolder, files(k).name));
        end
    elseif files(k).isdir && ~ismember(files(k).name, {'.', '..'})
        rmdir(fullfile(workspaceFolder, files(k).name), 's');
    end
end
