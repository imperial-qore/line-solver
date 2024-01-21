function jmtPath = jmtGetPath
% JMTPATH = JMTGETPATH

% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.

jmtPath = fileparts(which('JMT.jar'));
if isempty(jmtPath)
    line_printf('\nJava Modelling Tools cannot be found. LINE will try to download the latest JMT version (download approx. 30MB).\n')
    jmtSolverFolder = fileparts(mfilename('fullpath'));
%    if isdeployed
        m='Y';
%    else
%        m = input('Do you want to continue (download approx. 30MB), Y/N: ','s');
%    end
    if m=='Y'
        try
            line_printf('\nDownload started, please wait - this may take several minutes.')
            if exist('websave')==2
                if ispc
                    outfilename = websave([jmtSolverFolder,'\JMT.jar'],'http://jmt.sourceforge.net/latest/JMT.jar');
                else
                    outfilename = websave([jmtSolverFolder,'/JMT.jar'],'http://jmt.sourceforge.net/latest/JMT.jar');
                end
                line_printf('\nDownload completed. JMT jar now located at: %s',outfilename);
            else
                line_error(mfilename,'The MATLAB version is too old and JMT cannot be downloaded automatically. Please download http://jmt.sourceforge.net/latest/JMT.jar and put it under "bin\solvers\JMT\".\n');
            end
            jmtPath = fileparts(which('JMT.jar'));
        catch
            delete(which('JMT.jar'));
        end
    else
        if ispc
            line_error(mfilename,'Java Modelling Tools was not found. Please download http://jmt.sourceforge.net/latest/JMT.jar and put it under "bin\line\src\solvers\JMT\".\n');
        else
            line_error(mfilename,'Java Modelling Tools was not found. Please download http://jmt.sourceforge.net/latest/JMT.jar and put it under "bin/line/src/solvers/JMT/".\n');
        end
    end
end
end
