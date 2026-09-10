function savedfname = plot(model)
%PLOT Export a plot of an LQN model using lqn2ps
lastwd = pwd();
workspaceDir = fullfile(lineRootFolder, 'workspace', 'lqns');
cd(workspaceDir);

% Generate temp name with 'tmp_' prefix
rawname = tempname(workspaceDir);
[parentdir, filename] = fileparts(rawname);
fname = fullfile(parentdir, ['tmp_', filename]);
xmlFile = [fname, '.lqnx'];
model.writeXML(xmlFile, false);

system(sprintf('lqn2ps %s', xmlFile));
line_printf('Postscript file saved in: %s\n', [fname, '.ps']);

cd(lastwd);
savedfname = [fname, '.ps'];
end