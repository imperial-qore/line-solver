function solvePMIF(PMIFfilepath, options)
% SOLVEPMIF(PMIFFILEPATH, OPTIONS)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.

if nargin < 2
    options = [];
end
%% read options - set defaults if no options provided
[~, RT, RTrange, iter_max, solver, verbose] = parseOptions(options);

%% if this is a single file
if exist(PMIFfilepath, 'file') == 2 % it is a file
    allFilenames = {PMIFfilepath};
    [~, name, ext] = fileparts(PMIFfilepath);
    shortNames = {[name, ext]};
elseif exist(PMIFfilepath, 'dir') % it is a directory
    folderContents = dir(PMIFfilepath);
    names = {folderContents.name}';
    allFilenames = cell(0);
    shortNames = cell(0);
    for j = 1:size(names,1)
        if folderContents(j).isdir == 0 && strcmp(folderContents(j).name(end-3:end), '.xml')
            allFilenames{end+1,1} = fullfile(PMIFfilepath, names{j});
            shortNames{end+1,1} = names{j};
        end
    end
else
    line_printf(['Error: Input file path ', PMIFfilepath,' not found']);
    return;
end

for j = 1:size(allFilenames,1)
    myPath = allFilenames{j};
    % obtain CQN model from PMIF model
    sn = PMIF2QN(myPath, verbose);
    
    if ~isempty(sn)
        % check solver
        myRT=RT;
        % compute performance measures
        %iter_max = 1000;
        options=Solver.defaultOptions;
        [~, U, R, X, ~, RT_CDF, ~] = QN_fluid_analyzer(sn, [], [], [], myRT, RTrange, options);
        
        for i = 1:sn.nstations
            if sn.schedid(i)==SchedStrategy.ID_INF
                meanRT = sum(R([1:i-1 i+1:sn.nstations],:),1);
                break;
            end
        end
        
        % write results to output file
        outPath = [options.outputFolder, '/', shortNames{j}];
        writeXMLresults(outPath, '', sn, U, X, meanRT, R, [], RT_CDF, [], verbose );
    end
end
