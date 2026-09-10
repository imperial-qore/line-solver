function [kind, family, submethod] = resolveMethodToken(self, model, token)
% [KIND, FAMILY, SUBMETHOD] = RESOLVEMETHODTOKEN(MODEL, TOKEN)
%
% Split a method name into a selection intent, or into a method family and
% the submethod to be handed to it. A qualified method name 'family.submethod' keeps
% its submethod: dropping it would silently downgrade a pinned method to the
% family default.

if isempty(token)
    token = 'default';
end
token = char(token);
if strcmpi(token,'auto')
    token = 'default';
end

intents = SolverAUTO.selectionIntents();
if any(strcmpi(token, intents))
    kind = 'intent';
    family = lower(token);
    submethod = 'default';
    return
end

dotpos = find(token=='.', 1);
if isempty(dotpos)
    head = token;
    rest = '';
else
    head = token(1:dotpos-1);
    rest = token(dotpos+1:end);
end

fam = SolverAUTO.familyAlias(head);
if ~isempty(fam)
    kind = 'family';
    family = fam;
    if isempty(rest)
        % A bare family name means its default method; for bounds the default
        % is the composite tightest-of-all family rather than a single one.
        if strcmp(fam,'ba')
            submethod = 'auto';
        else
            submethod = 'default';
        end
    else
        submethod = rest;
    end
    return
end

% Unqualified algorithm name, e.g. 'comom' or 'gb.upper'. The family that
% declares it owns it, so the name table stays in the families themselves.
fam = self.familyDeclaringMethod(model, token);
if ~isempty(fam)
    kind = 'family';
    family = fam;
    submethod = token;
    return
end

line_error(mfilename, sprintf(['Unrecognized method ''%s''. Valid tokens are a selection intent (%s), ', ...
    'a method family (%s), or a qualified method name such as ''nc.comom''.'], ...
    token, strjoin(intents, ', '), strjoin(SolverAUTO.familyNames(), ', ')));
end
