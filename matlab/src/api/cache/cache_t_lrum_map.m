%{ @file cache_t_lrum_map.m
 %  @brief Characteristic times for LRU(m) caches under MAP request streams
 %
 %  @author LINE Development Team
%}

%{
 % @brief Characteristic times for the LRU(m)-MAP TTL approximation
 %
 % @details
 % TTL approximation of LRU(m) with per-item Markovian arrival processes
 % (Gast and Van Houdt, Performance Evaluation 2017, Section 3.1.2). Each
 % item is modeled by an embedded Markov chain over (list, phase) states;
 % the level vectors obey pi_l = pi_0 prod_s R_s with the R-recursions of
 % eqs. (6)-(7) and pi_0 the left Perron vector of R_1 expm(D0 T_1). The
 % characteristic times T_1..T_h equate the expected occupancy of each
 % list to its capacity; the joint root is found by fsolve on log(T).
 %
 % @par Syntax:
 % @code
 % t = cache_t_lrum_map(D0c, D1c, m)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>D0c<td>Cell (1,n); D0c{k} is the (d,d) hidden-transition matrix of item k
 % <tr><td>D1c<td>Cell (1,n); D1c{k} is the (d,d) arrival matrix of item k
 % <tr><td>m<td>Cache capacity vector (1,h)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>t<td>Characteristic time of each list (1,h)
 % </table>
%}
function t = cache_t_lrum_map(D0c, D1c, m)
n = numel(D0c);
h = numel(m);

    function F = capres(y)
        T = exp(y);
        occ = zeros(1, h);
        for k = 1:n
            [~, occk] = cache_lrum_map_levelstats(D0c{k}, D1c{k}, T);
            occ = occ + occk;
        end
        F = occ - m(:)';
    end

options = optimoptions('fsolve', 'MaxIter', 1e4, 'MaxFunEvals', 1e5, 'Display', 'off');
y = fsolve(@capres, zeros(1, h), options);
t = exp(y);
end
