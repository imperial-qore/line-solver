%{ @file ba_resolve_method.m
 %  @brief Resolves the SolverBA method aliases to a concrete bound name
 %
 %  @author LINE Development Team
%}

%{
 % @brief Resolves 'default', 'auto', 'qr' and 'lr' to the method they name
 %
 % @details
 % 'default' selects the tightest noniterative upper bound (geometric); bare
 % 'auto' is the AUTO composite's upper side; 'qr' is the friendly alias of
 % the QRF quadratic reduction bound; bare 'lr' is the LP-based Linear
 % Reduction upper bound (Mapqn_bnd_lr_pf, simplex), NOT an alias of
 % 'qrf.mmi.linear', which is "linear" only in its explicit Aeq/beq
 % constraint representation while its objective is the nonlinear MMI
 % mutual information.
 %
 % Any other name is returned unchanged. Kept as one function so that
 % runAnalyzer, listValidMethods and getBounds all judge a method by the same
 % resolved name: gating on the raw 'default' would let the geometric bound
 % through every check written against 'gb.upper'.
 %
 % BARE 'auto' IS DELIBERATELY ABSENT FROM SolverBA.listAllMethods, and the
 % arm here is live rather than dead. SolverBA is the one Network solver that
 % does not pass through NetworkSolver.runAnalyzerChecks, so no name gate ever
 % reads a list before this resolution runs: runAnalyzer dispatches on the
 % RESOLVED 'auto.upper', which listAllMethods does carry. Three production
 % callers set it -- SolverAUTO's 'bound' intent (SolverAUTO.m), its baSolver
 % fallback, and resolveMethodToken's bare 'ba' token, all added with 'bound'
 % in 00fc0acaf. Listing it would put it in listValidMethods too (that list is
 % a narrowing of listAllMethods), which enumerates into the qn/jqn sanity
 % baselines; the alias needs no listing to work, so it stays unlisted.
 %
 % @par Syntax:
 % @code
 % method = ba_resolve_method(method)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>method<td>SolverBA bound method name as given
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>method<td>The resolved method name
 % </table>
%}
function method = ba_resolve_method(method)
if strcmp(method,'default')
    method = 'gb.upper';
elseif strcmp(method,'auto')
    method = 'auto.upper';
elseif strcmp(method,'lr')
    method = 'lr.upper';
elseif strcmp(method,'qr')
    method = 'qrf.mmi';
end
end
