function model = LINE2JAVA(model, filename)
% MODEL = LINE2JAVA(MODEL, FILENAME)

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.
if nargin>=2 %exist('filename','var')
    fid = fopen(filename,'w'); % discard
    if isa(model,'Network')
        QN2JAVA(model, model.getName(), fid);
    elseif isa(model,'LayeredNetwork')
        LQN2JAVA(model, model.getName(), fid);
    end
    fclose(fid);
else
    if isa(model,'Network')
        QN2JAVA(model, model.getName(), 1);
    elseif isa(model,'LayeredNetwork')
        LQN2JAVA(model, model.getName(), 1);
    end
end
end
