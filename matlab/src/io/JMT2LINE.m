function model = JMT2LINE(filename,modelName)
% MODEL = JMT2LINE(FILENAME,MODELNAME)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[~,~,fext] = fileparts(filename);
switch fext
    case '.jmva'
        line_error(mfilename,'JMVA files not yet supported.');
        if nargin<2
            model = JMVA2LINE(filename);
        else
            model = JMVA2LINE(filename, modelName);
        end
    case {'.jsim','.jsimg','.jsimw'}
        % create network
        if nargin<2
            model = JSIM2LINE(filename);
        else
            model = JSIM2LINE(filename, modelName);
        end
    otherwise
        line_warning(mfilename,'The file has unknown extension, trying to parse as a JSIMG file.\n');
        if nargin<2
            model = JSIM2LINE(filename);
        else
            model = JSIM2LINE(filename, modelName);
        end
end

end
