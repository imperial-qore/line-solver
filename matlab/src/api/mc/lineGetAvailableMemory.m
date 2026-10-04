function bytes = lineGetAvailableMemory()
% LINEGETAVAILABLEMEMORY  Portable available-physical-memory probe (bytes).
%
% Routes through the JVM OperatingSystemMXBean that MATLAB ships on every
% platform, so the same code path returns a valid figure on Windows, macOS
% and Linux without parsing /proc. Falls back to a conservative constant if
% the bean or method is unavailable (e.g. MATLAB started with -nojvm).
%
% Mirrors jline.solvers.ctmc.MemoryGuard.getAvailableMemoryBytes (JAR) and
% line_solver ... memory_guard.get_available_memory_bytes (Python native).

FALLBACK_BYTES = 1 * 1024^3;  % 1 GB
bytes = FALLBACK_BYTES;

try
    osb = java.lang.management.ManagementFactory.getOperatingSystemMXBean();
    v = -1;
    % JDK <= 13: getFreePhysicalMemorySize; JDK >= 14: getFreeMemorySize
    try
        v = osb.getFreePhysicalMemorySize();
    catch
        try
            v = osb.getFreeMemorySize();
        catch
            v = -1;
        end
    end
    if ~isempty(v) && v > 0
        bytes = double(v);
    end
catch
    % leave fallback
end
end
