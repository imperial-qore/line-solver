function LINE2MATLAB(model, filename)
% MODEL = LINE2MATLAB(MODEL, FILENAME)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if ~isa(model,'MNetwork') && ~isa(model,'LayeredNetwork')
    line_error(mfilename, sprintf('LINE2MATLAB supports a Network or a LayeredNetwork, got a %s.', class(model)));
end
if nargin>=2 %exist('filename','var')
    fid = fopen(filename,'w');
    if isa(model,'MNetwork')
        QN2MATLAB(model, model.getName(), fid);
    elseif isa(model,'LayeredNetwork')
        LQN2MATLAB(model, model.getName(), fid);
    end
    fclose(fid);
else
    if isa(model,'MNetwork')
        QN2MATLAB(model, model.getName(), 1);
    elseif isa(model,'LayeredNetwork')
        LQN2MATLAB(model, model.getName(), 1);
    end
end
end
