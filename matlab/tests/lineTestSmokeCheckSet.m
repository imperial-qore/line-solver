function lineTestSmokeCheckSet(testName, solverNames, fieldNames, valuesByField)
% LINETESTSMOKECHECKSET  Assert that a dispatch run returned usable numbers.
%
% The completion half of smoke mode (see LINETESTSMOKEMODE). It is handed, per
% metric, the value each solver produced -- valuesByField{f}{s} for field
% fieldNames{f} and solver solverNames{s} -- and fails the test unless every one
% of them is a real, finite, non-negative numeric array.
%
% AN EMPTY ENTRY IS NOT A FAILURE. A solver that refused the model by featset is
% caught by the caller and leaves its slot empty on purpose, which is the one way
% a dispatch backend is allowed to answer with nothing. A run in which EVERY slot
% is empty is a failure, though: that is a suite that tested nothing while
% reporting success, the exact silence smoke mode exists to break.
%
% No upper bound is asserted on utilization: a multiserver station is measured
% per station, not per server, so a legitimate value there exceeds one.

if ~iscell(fieldNames) || ~iscell(valuesByField) || numel(fieldNames) ~= numel(valuesByField)
    error('lineTest:smoke', '%s: %d field names against %d value sets.', ...
        testName, numel(fieldNames), numel(valuesByField));
end

anyValue = false;

for f = 1:numel(fieldNames)
    field = fieldNames{f};
    values = valuesByField{f};
    if ~iscell(values)
        values = {values};
    end

    for s = 1:numel(values)
        val = values{s};
        if s <= numel(solverNames)
            who = char(solverNames{s});
        else
            who = sprintf('solver %d', s);
        end

        if isempty(val)
            continue    % featset refusal: the caller has already allowed it
        end
        anyValue = true;

        if ~isnumeric(val) && ~islogical(val)
            error('lineTest:smoke', '%s: %s returned a %s for %s, not numbers.', ...
                testName, who, class(val), field);
        end
        val = double(val(:));
        if ~isreal(val)
            error('lineTest:smoke', '%s: %s returned a complex %s.', testName, who, field);
        end
        bad = find(~isfinite(val), 1);
        if ~isempty(bad)
            error('lineTest:smoke', '%s: %s has a non-finite %s at entry %d (%g).', ...
                testName, who, field, bad, val(bad));
        end
        % Queue lengths, utilizations, response times and throughputs are all
        % non-negative; the tolerance is there for a solver that lands on zero
        % from below rather than to accept a negative answer.
        bad = find(val < -1e-9, 1);
        if ~isempty(bad)
            error('lineTest:smoke', '%s: %s has a negative %s at entry %d (%g).', ...
                testName, who, field, bad, val(bad));
        end
    end
end

if ~anyValue
    error('lineTest:smoke', '%s: every solver returned nothing, so the test asserted nothing.', ...
        testName);
end
end
