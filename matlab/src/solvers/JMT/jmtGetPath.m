function jmtPath = jmtGetPath
% JMTPATH = JMTGETPATH

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% Get path to common folder
% This function is at dev/src/solvers/JMT/jmtGetPath.m
% Navigate up to the line-dev.git root directory
current_dir = fileparts(mfilename('fullpath'));           % dev/src/solvers/JMT
root_dir = fileparts(fileparts(fileparts(fileparts(current_dir)))); % line-dev.git
common_dir = fullfile(root_dir, 'common');
jmt_jar = fullfile(common_dir, 'JMT.jar');

if exist(jmt_jar, 'file')
    jmtPath = common_dir;
else
    line_printf('\nJava Modelling Tools cannot be found. LINE will try to download the latest JMT version (download approx. 50MB).\n')
    m='Y';
    if m=='Y'
        try
            line_printf('\nDownload started, please wait - this may take several minutes.')
            if exist('websave')==2
                % Check if file exists but is locked/permission issue
                if exist(jmt_jar, 'file')
                    fprintf('JMT.jar already exists at %s, attempting to use existing file.\n', jmt_jar);
                    jmtPath = common_dir;
                    return;
                end
                
                % Ensure common directory exists
                if ~exist(common_dir, 'dir')
                    mkdir(common_dir);
                end
                
                jmt_url = 'https://line-solver.sourceforge.net/latest/JMT.jar';
                outfilename = websave(jmt_jar, jmt_url);
                line_printf('\nDownload completed. JMT jar now located at: %s',outfilename);
            else
                line_error(mfilename,'The MATLAB version is too old and JMT cannot be downloaded automatically. Please download https://line-solver.sourceforge.net/latest/JMT.jar and put it in the common folder.');
            end
            jmtPath = common_dir;
        catch ME
            if exist(jmt_jar, 'file')
                delete(jmt_jar);
            end
            error_msg = sprintf('Failed to download JMT: %s\nPlease manually download https://line-solver.sourceforge.net/latest/JMT.jar and place it in %s', ME.message, common_dir);
            line_error(mfilename, error_msg);
        end
    else
        error_msg = sprintf('Java Modelling Tools was not found. Please download https://line-solver.sourceforge.net/latest/JMT.jar and put it in %s\n', common_dir);
        line_error(mfilename, error_msg);
    end
end
end
