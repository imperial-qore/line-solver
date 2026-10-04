function javaExe = line_java_exe()
% JAVAEXE = LINE_JAVA_EXE()
%
% Resolve the Java launcher used to spawn JMT (and any other JVM tool driven
% out of MATLAB). Returns '' when no JVM is reachable, so a caller can raise a
% diagnosis of its own instead of letting the failure surface as a raw
% java.io.IOException ("Cannot run program ""java"": CreateProcess error=2").
%
% Search order:
%   1. LINE_JAVA         - full path of a launcher chosen by the user
%   2. JAVA_HOME/bin     - the conventional JDK/JRE variable
%   3. MATLAB_JAVA/bin   - the JVM MATLAB itself was pointed at
%   4. the JRE bundled with MATLAB (layout varies across releases)
%   5. "java" on PATH
%
% Step 4 is what keeps the JMT viewers usable on a host with no system-wide
% Java, the common case on Windows: MATLAB ships its own JRE, so a missing
% PATH entry is not a missing JVM.
%
% The result is cached per session: the PATH probe costs a process launch, and
% a JVM neither appears nor vanishes mid-session.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

persistent cached
if ~isempty(cached)
    javaExe = cached{1};
    return
end

if ispc
    exeName = 'java.exe';
else
    exeName = 'java';
end

cands = {};
envJava = getenv('LINE_JAVA');
if ~isempty(envJava)
    cands{end+1} = envJava;
end
javaHome = getenv('JAVA_HOME');
if ~isempty(javaHome)
    cands{end+1} = fullfile(javaHome, 'bin', exeName);
end
matlabJava = getenv('MATLAB_JAVA');
if ~isempty(matlabJava)
    cands{end+1} = fullfile(matlabJava, 'bin', exeName);
end
cands{end+1} = fullfile(matlabroot, 'sys', 'java', 'jre', computer('arch'), 'jre', 'bin', exeName);
cands = [cands, lineJreGlob(exeName)];

javaExe = '';
for i = 1:numel(cands)
    if exist(cands{i}, 'file') == 2
        javaExe = cands{i};
        cached = {javaExe};
        return
    end
end

% Last resort: whatever PATH resolves. A launcher that is present but broken
% (wrong architecture, truncated install) fails here too, which is the point.
[st, ~] = system([exeName, ' -version 2>&1']);
if st == 0
    javaExe = exeName;
end
cached = {javaExe};
end


function paths = lineJreGlob(exeName)
% PATHS = LINEJREGLOB(EXENAME)
% Bundled-JRE candidates found by walking matlabroot rather than by assuming
% one layout: the directory moved between releases and is absent in some.
paths = {};
jreRoot = fullfile(matlabroot, 'sys', 'java', 'jre');
if exist(jreRoot, 'dir') ~= 7
    return
end
entries = dir(jreRoot);
for i = 1:numel(entries)
    if ~entries(i).isdir || any(strcmp(entries(i).name, {'.', '..'}))
        continue
    end
    archDir = fullfile(jreRoot, entries(i).name);
    paths{end+1} = fullfile(archDir, 'jre', 'bin', exeName); %#ok<AGROW>
    paths{end+1} = fullfile(archDir, 'bin', exeName); %#ok<AGROW>
end
end
