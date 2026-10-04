function LINE2JAVA(model, filename)
% LINE2JAVA(MODEL, FILENAME)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if ~isa(model,'MNetwork') && ~isa(model,'LayeredNetwork')
    line_error(mfilename, sprintf('LINE2JAVA supports a Network or a LayeredNetwork, got a %s.', class(model)));
end
if nargin>=2 %exist('filename','var')
    fid = fopen(filename,'w');
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
