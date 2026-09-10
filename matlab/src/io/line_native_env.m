function prefix = line_native_env()
% PREFIX = LINE_NATIVE_ENV()
%
% Shell prefix that restores a SYSTEM loader environment for a native binary
% spawned out of MATLAB (line-cli, the ldes native binary).
%
% MATLAB exports its own bin/glnxa64 and sys/os/glnxa64 to every child process,
% and the libstdc++ shipped there is older than the one a current toolchain
% links against. A binary this checkout compiled then fails to load:
%
%   common/line-cli: .../sys/os/glnxa64/libstdc++.so.6: version
%   `GLIBCXX_3.4.32' not found (required by common/line-cli)
%
% which surfaces as "line-cli exited with code 1" and, one level up, as every
% lang='cpp' row of the model reporting "solver X missing from output". Neither
% the binary nor the model is at fault, so the fix belongs here.
%
% Only the entries under matlabroot (and MATLAB's per-user .MathWorks tree) are
% dropped: a path the user set before starting MATLAB is theirs and is kept.
% A JVM runner must NOT be wrapped -- it is MATLAB's own JRE and wants those
% directories.
%
% Returns '' on Windows and macOS, where no such shadowing occurs.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

prefix = '';
if ~isunix || ismac
    return
end

ldPath = getenv('LD_LIBRARY_PATH');
if isempty(ldPath)
    return
end

mlRoot = matlabroot;
entries = strsplit(ldPath, ':');
keep = {};
for i = 1:numel(entries)
    e = entries{i};
    if isempty(e) || strncmp(e, mlRoot, numel(mlRoot)) || contains(e, '/.MathWorks/')
        continue
    end
    keep{end+1} = e; %#ok<AGROW>
end
% An emptied variable is still exported as empty: omitting it would let the
% child inherit MATLAB's value instead of the stripped one.
prefix = sprintf('env LD_LIBRARY_PATH="%s" ', strjoin(keep, ':'));
end
