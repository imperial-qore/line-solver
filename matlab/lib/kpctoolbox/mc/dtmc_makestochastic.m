%{ @file dtmc_makestochastic.m
 %  @brief Normalizes a matrix to be a valid stochastic matrix
 %
 %  @author LINE Development Team
%}

%{
 % @brief Normalizes a matrix to be a valid stochastic matrix
 %
 % @details
 % Rescales rows to sum to 1. Rows with zero sum are set to have 1 at the diagonal.
 %
 % @par Syntax:
 % @code
 % P = dtmc_makestochastic(P)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>P<td>Input matrix
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>P<td>Valid stochastic transition matrix
 % </table>
%}
function P=dtmc_makestochastic(P)
for i=1:size(P,1)
    if sum(P(i,:))>0
        P(i,:)=P(i,:)/sum(P(i,:));
        P(i,i)=min(max(0,1-(sum(P(i,:))-P(i,i))),1);
    else
        P(i,:)=0;
        P(i,i)=1;
    end
end
end