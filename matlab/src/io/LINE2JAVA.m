function LINE2JAVA(model, filename)
% LINE2JAVA(MODEL, FILENAME)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin>=2 %exist('filename','var')
    fid = fopen(filename,'w'); % discard
    if isa(model,'MNetwork')
        QN2JAVA(model, model.getName(), fid);
    elseif isa(model,'LayeredNetwork')
        LQN2JAVA(model, model.getName(), fid);
    end
    fclose(fid);
else
    if isa(model,'MNetwork')
        QN2JAVA(model, model.getName(), 1);
    elseif isa(model,'LayeredNetwork')
        LQN2JAVA(model, model.getName(), 1);
    end
end
end
