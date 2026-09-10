function [tf, freeBytes, requiredBytes] = lineDockerHasStorageFor(image)
% LINEDOCKERHASSTORAGEFOR True if the Docker storage location has room for IMAGE.
%
%   [TF, FREEBYTES, REQUIREDBYTES] = LINEDOCKERHASSTORAGEFOR(IMAGE)
%
% Compares the free space on the filesystem backing the Docker root dir against
% max(2 GiB floor, estimated on-disk image size). Returns TF=true when free
% space cannot be determined (do not block the pull). The required floor can be
% overridden with the LINE_DOCKER_MIN_FREE_BYTES environment variable.
%
% Shared by the LQNS/QNS resolver (lineDockerImage) and the JMT wrapper so the
% pre-pull storage guard has a single implementation. Mirrors the Java
% DockerImage.hasStorageFor and the Python docker_util.has_storage_for.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

freeBytes = localFreeBytes();
requiredBytes = localRequiredBytes(image);
if freeBytes < 0
    tf = true; % could not determine; do not block the pull
else
    tf = (freeBytes >= requiredBytes);
end
end

function bytes = localFreeBytes()
% Usable bytes on the filesystem backing the Docker root dir; -1 if unknown.
bytes = -1;
[st, out] = unix('docker info --format "{{.DockerRootDir}}" 2>/dev/null');
root = '';
if st == 0
    root = strtrim(out);
end
if isempty(root)
    root = '/var/lib/docker';
end
% The root dir may not exist for / be readable by this user: walk up to an
% existing ancestor so getUsableSpace() reports a real filesystem.
while ~isempty(root) && exist(root, 'dir') ~= 7
    parent = fileparts(root);
    if strcmp(parent, root)
        break;
    end
    root = parent;
end
if isempty(root) || exist(root, 'dir') ~= 7
    root = '/';
end
try
    f = java.io.File(root);
    usable = f.getUsableSpace();
    if usable > 0
        bytes = double(usable);
    end
catch
    bytes = -1;
end
end

function bytes = localRequiredBytes(image)
DEFAULT_MIN = 2 * 1024 * 1024 * 1024; % 2 GiB floor
ov = getenv('LINE_DOCKER_MIN_FREE_BYTES');
if ~isempty(ov)
    v = str2double(ov);
    if ~isnan(v) && v > 0
        bytes = v;
        return;
    end
end
bytes = max(DEFAULT_MIN, localEstimateBytes(image));
end

function bytes = localEstimateBytes(image)
% On-disk estimate (compressed layer sizes x3), or 0 if it cannot be determined.
bytes = 0;
[st, out] = unix(['docker manifest inspect ', image, ' 2>/dev/null']);
if st ~= 0 || isempty(out)
    return;
end
tok = regexp(out, '"size"\s*:\s*(\d+)', 'tokens');
total = 0;
for i = 1:numel(tok)
    total = total + str2double(tok{i}{1});
end
if total > 0
    bytes = total * 3;
end
end
