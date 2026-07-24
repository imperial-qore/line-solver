function build_lsoda_mex()
% BUILD_LSODA_MEX  Compile the LSODA MEX solver
%
%   build_lsoda_mex()
%
%   Compiles liblsoda C sources into a MEX function lsoda_mex.
%   The compiled MEX file is placed in the same directory as this script.
%
%   Requirements: a C compiler configured via mex -setup

thisDir = fileparts(mfilename('fullpath'));
srcDir  = fullfile(thisDir, 'src');

% List all C source files (excluding ewset.c and printcf.c)
srcFiles = {
    fullfile(srcDir, 'lsoda_mex.c')
    fullfile(srcDir, 'lsoda.c')
    fullfile(srcDir, 'stoda.c')
    fullfile(srcDir, 'correction.c')
    fullfile(srcDir, 'prja.c')
    fullfile(srcDir, 'solsy.c')
    fullfile(srcDir, 'intdy.c')
    fullfile(srcDir, 'cfode.c')
    fullfile(srcDir, 'methodswitch.c')
    fullfile(srcDir, 'orderswitch.c')
    fullfile(srcDir, 'corfailure.c')
    fullfile(srcDir, 'scaleh.c')
    fullfile(srcDir, 'common.c')
    fullfile(srcDir, 'vmnorm.c')
    fullfile(srcDir, 'fnorm.c')
    fullfile(srcDir, 'daxpy.c')
    fullfile(srcDir, 'ddot.c')
    fullfile(srcDir, 'dgefa.c')
    fullfile(srcDir, 'dgesl.c')
    fullfile(srcDir, 'dscal.c')
    fullfile(srcDir, 'idamax.c')
    fullfile(srcDir, 'strdup_printf.c')
};

fprintf('Building lsoda_mex...\n');

mex('-O', ...
    ['-I' srcDir], ...
    '-output', fullfile(thisDir, 'lsoda_mex'), ...
    srcFiles{:});

fprintf('lsoda_mex built successfully.\n');
end
