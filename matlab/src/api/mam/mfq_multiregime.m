%{ @file mfq_multiregime.m
 %  @brief Multi-regime feedback fluid queue solver
 %
 %  @author LINE Development Team
%}

%{
 % @brief Solves a multi-regime feedback Markovian fluid queue.
 %
 % @details
 % Thin LINE wrapper around the multiregime routine implementing the method of
 % H. E. Kankaya and N. Akar, "Solving Multi-Regime Feedback Fluid Queues".
 % The generator and drift rates are regime dependent, with separate boundary
 % behaviour (feedback) generators/rates at each threshold.
 %
 % @par Syntax:
 % @code
 % [pdf,pdfd,cdf,cdfm] = mfq_multiregime(Q,R,Qt,Rt,T,pdfpoints,cdfpoints)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>Cell of (N,N) generators, one per regime k=1..K
 % <tr><td>R<td>Cell of (N,N) diagonal drift-rate matrices per regime
 % <tr><td>Qt<td>Cell of boundary generators, k=0..K
 % <tr><td>Rt<td>Cell of boundary drift-rate matrices, k=0..K
 % <tr><td>T<td>Vector of regime thresholds (length K)
 % <tr><td>pdfpoints<td>Levels at which to evaluate density / density derivative
 % <tr><td>cdfpoints<td>Levels at which to evaluate the distribution function
 % </table>
 %
 % @par Returns:
 % pdf, pdfd (density derivative), cdf P(X<p), and cdfm P(X<=p).
%}
function [pdf,pdfd,cdf,cdfm] = mfq_multiregime(Q,R,Qt,Rt,T,pdfpoints,cdfpoints)
[pdf,pdfd,cdf,cdfm] = multiregime(Q,R,Qt,Rt,T,pdfpoints,cdfpoints);
end
