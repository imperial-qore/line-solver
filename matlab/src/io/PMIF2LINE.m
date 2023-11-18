function model = PMIF2LINE(filename,modelName)
% MODEL = PMIF2LINE(FILENAME,MODELNAME)

% Copyright (c) 2012-2023, Imperial College London
% All rights reserved.
verbose = false;
sn = PMIF2QN(filename,verbose);
model = QN2LINE(sn, modelName);
end
