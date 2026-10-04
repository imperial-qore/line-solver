%{ @file mfq_ld_solve.m
 %  @brief Solves first/second-order level-dependent (multi-regime) fluid queues
 %
 %  @author LINE Development Team
%}

%{
 % @brief Matrix-exponential solution of a multi-regime Markovian fluid queue.
 %
 % @details
 % Thin LINE wrapper around SecondOrderLevelDependentFluidSolve. The
 % generator, drift and (optionally) variance change at threshold fluid
 % levels, yielding a piecewise-homogeneous first- or second-order (Brownian)
 % fluid queue. Setting the variance cells S to zero reduces the model to
 % first order. The returned matrix-exponential building blocks are consumed
 % by mfq_ld_distr / mfq_ld_mean.
 %
 % @par Syntax:
 % @code
 % [masses,iniF,KF,cloF,iniB,KB,cloB] = mfq_ld_solve(Q,R,S,T)
 % [...] = mfq_ld_solve(Q,R,S,T,boundaryL,boundaryU,Qt,prec)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>Cell of (N,N) generators, one per regime
 % <tr><td>R<td>Cell of (N,N) diagonal drift-rate matrices per regime
 % <tr><td>S<td>Cell of (N,N) diagonal variance matrices per regime (0 = first order)
 % <tr><td>T<td>Vector of regime thresholds (length K)
 % <tr><td>boundaryL/U<td>(Optional) per-state boundary flag (0 reflective, 1 absorbing)
 % <tr><td>Qt<td>(Optional) cell of boundary generators
 % <tr><td>prec<td>(Optional) numerical precision (default 1e-14)
 % </table>
 %
 % @par Returns:
 % masses (K+1 point-mass vectors) and per-regime forward/backward
 % matrix-exponential parameters iniF,KF,cloF,iniB,KB,cloB.
%}
function varargout = mfq_ld_solve(varargin)
[varargout{1:nargout}] = SecondOrderLevelDependentFluidSolve(varargin{:});
end
