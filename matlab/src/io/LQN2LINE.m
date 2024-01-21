function model = LQN2LINE(filename, modelName)
% MODEL = LQN2LINE(FILENAME, MODELNAME)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

[~,~,fext] = fileparts(filename);
switch fext
    case {'lqnx','xml','.xml','.lqnx'}
        % create network        
        model = LayeredNetwork.parseXML(filename, false);
        model.setName(modelName);
    otherwise
        line_error(mfilename, 'File extension must be .lqnx or .xml');
end

end
