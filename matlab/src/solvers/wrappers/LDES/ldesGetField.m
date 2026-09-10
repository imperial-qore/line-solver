function v = ldesGetField(s, f, dflt)
% V = LDESGETFIELD(S, F, DFLT)
% Value of field F in struct S, or DFLT if S is not a struct or lacks F.
% Shared helper for the fully JSON-mediated LDES wrapper methods.
if isstruct(s) && isfield(s, f)
    v = s.(f);
else
    v = dflt;
end
end
