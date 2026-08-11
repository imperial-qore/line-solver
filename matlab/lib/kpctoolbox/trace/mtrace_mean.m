function m = mtrace_mean(trace, ntypes, type)
% [m]=mtrace_mean(trace, ntypes, type)
%
% DESCRIPTION
% Compute the mean of a trace, divided by types
%
% PARAMETERS
% trace  - the array containing the trace data
% ntypes - the number of different types
% type   - an array indicating the type of each element in the trace (0-indexed)
%
% RETURN
% m - a column vector containing the mean values for each type

% Initialize mean vector
m = zeros(ntypes, 1);

% Compute mean for each type
for c = 0:(ntypes-1)
    % Find indices where type equals c
    idx = (type == c);
    if any(idx)
        % Compute mean only for elements of this type
        m(c+1) = mean(trace(idx));
    else
        % If no elements of this type, set to NaN
        m(c+1) = NaN;
    end
end

end