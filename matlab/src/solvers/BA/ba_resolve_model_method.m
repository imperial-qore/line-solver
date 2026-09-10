%{ @file ba_resolve_model_method.m
 %  @brief Resolves a SolverBA method name to the method it runs as on a model
 %
 %  @author LINE Development Team
%}

%{
 % @brief The concrete bound method NAME runs as on the model SN
 %
 % @details
 % Two resolutions in one place, so that runAnalyzer, SolverBA.resolveMethod,
 % BA_METHOD_REFUSAL and SolverBA.getMethodFeatureSet all judge a name by the
 % same target. First the model-free aliases of BA_RESOLVE_METHOD ('default'
 % -> gb.upper, 'auto' -> auto.upper, 'qr' -> qrf.mmi, 'lr' -> lr.upper). Then
 % the finite-buffer routing of BA_BLOCKING_DEFAULT: on a model whose buffer
 % BINDS, 'default', 'auto' and 'auto.upper' run as 'qrf.bas', the one upper
 % bound that models the buffer, when the model has the shape that bound needs.
 %
 % WHY is the reason the routing did NOT apply -- the delay station, the
 % multiserver station, the second binding buffer or the oversized enumeration
 % that BA_BLOCKING_DEFAULT names -- and is empty when it applied, when the
 % name is not one of the three routable ones, or when the model is unblocked.
 % BA_METHOD_REFUSAL appends it to the blocking refusal so that a caller is
 % told what about THIS model rules the blocking bound out.
 %
 % Only the UPPER side routes: the analyzer solves 'qrf.bas' in the 'max'
 % direction alone, so 'auto.lower' keeps refusing rather than being answered
 % with the wrong side.
 %
 % @par Syntax:
 % @code
 % [method, why] = ba_resolve_model_method(sn, method)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>sn<td>Network structure
 % <tr><td>method<td>SolverBA bound method name as given
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>method<td>The method the name runs as on this model
 % <tr><td>why<td>Reason the finite-buffer routing did not apply, '' otherwise
 % </table>
%}
function [method, why] = ba_resolve_model_method(sn, method)
why = '';
routable = any(strcmp(method, {'default','auto','auto.upper'}));
method = ba_resolve_method(method);
if routable && ba_ignores_blocking(method) && sn_has_blocking(sn)
    [alt, why] = ba_blocking_default(sn);
    if ~isempty(alt)
        method = alt;
    end
end
end
