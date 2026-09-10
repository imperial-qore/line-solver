%{ @file mfq_prio_queue.m
 %  @brief Performance measures of a continuous-time fluid priority queue
 %
 %  @author LINE Development Team
%}

%{
 % @brief Analyzes an MMAP[K]/PH[K]/1-type fluid priority queue.
 %
 % @details
 % Thin LINE wrapper around the BUTools-family routine FluidPrioQueue,
 % implementing the method of G. Horvath, "Efficient analysis of the
 % MMAP[K]/PH[K]/1 priority queue", EJOR 246(1):128-139, 2015. A background
 % Markov chain with generator Q modulates the per-class fluid input rates
 % (matrix R, one row per priority class) and the fluid is drained at the
 % constant service rate d, higher-priority fluid first.
 %
 % @par Syntax:
 % @code
 % varargout = mfq_prio_queue(Q, R, d, ...)
 % [flM] = mfq_prio_queue(Q, R, d, 'flMoms', n)
 % [cdf] = mfq_prio_queue(Q, R, d, 'stDistr', points, 'erlMaxOrder', L)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>(N,N) generator of the modulating Markov chain
 % <tr><td>R<td>(K,N) per-class fluid input rates in the background states
 % <tr><td>d<td>Constant fluid service rate (scalar, positive)
 % <tr><td>...<td>Measure/option pairs: 'flMoms','flDistr','stMoms','stDistr',
 %               'prec','erlMaxOrder','classes' (see FluidPrioQueue)
 % </table>
 %
 % @par Returns:
 % One output per requested performance measure; each column corresponds to a
 % priority class.
%}
function varargout = mfq_prio_queue(Q, R, d, varargin)
[varargout{1:nargout}] = FluidPrioQueue(Q, R, d, varargin{:});
end
