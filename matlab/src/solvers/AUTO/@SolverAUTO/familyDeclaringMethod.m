function family = familyDeclaringMethod(self, model, methodName)
% FAMILY = FAMILYDECLARINGMETHOD(MODEL, METHOD NAME)
%
% Find the method family that declares an unqualified algorithm name, by
% asking each family for its own method list. Keeping the question there
% avoids a second copy of the name table in SolverAUTO, which would drift.

family = '';
order = SolverAUTO.familyNames();
probeOptions = Solver.defaultOptions;
probeOptions.verbose = 0;
for f = 1:length(order)
    try
        probe = self.buildFamilySolver(order{f}, model, probeOptions);
        declared = probe.listValidMethods();
    catch
        % A family that cannot even be instantiated on this model cannot own
        % the method name; the next one is asked instead.
        continue
    end
    if any(strcmpi(methodName, declared))
        family = order{f};
        return
    end
end
end
