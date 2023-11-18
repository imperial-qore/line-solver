function model = LINE2SCRIPT(model, filename)
% MODEL = LINE2SCRIPT(MODEL, FILENAME)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.
if nargin>=2 %exist('filename','var')
    fid = fopen(filename,'w'); % discard
    if isa(model,'Network')
        QN2SCRIPT(model, model.getName(), fid);
    elseif isa(model,'LayeredNetwork')
        LQN2SCRIPT(model, model.getName(), fid);
    end
    fclose(fid);
else
    if isa(model,'Network')
        QN2SCRIPT(model, model.getName(), 1);
    elseif isa(model,'LayeredNetwork')
        LQN2SCRIPT(model, model.getName(), 1);
    end
end
end
