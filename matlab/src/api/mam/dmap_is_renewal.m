%{ @file dmap_is_renewal.m
 %  @brief Tests whether a DMAP is a renewal process
 %
 %  @author LINE Development Team
%}

%{
 % @brief Tests whether a DMAP is a renewal process
 %
 % @details
 % A DMAP renews at every event exactly when D1 has rank one, i.e. D1 = a*alpha
 % with a the absorption vector of D0. Interevent times are then i.i.d. and the
 % process is the DMAP form of a discrete phase-type law, which is what lets
 % the caller use the cheaper Q_DT_PH_PH_1 entry point.
 %
 % @par Syntax:
 % @code
 % bool = dmap_is_renewal(DMAP)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>DMAP<td>Cell {D0, D1}
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>bool<td>True when D1 has rank one
 % </table>
%}
function bool = dmap_is_renewal(DMAP)

D1 = DMAP{2};
if size(D1,1) == 1
    bool = true;
    return;
end
bool = rank(D1, GlobalConstants.FineTol) <= 1;

end
