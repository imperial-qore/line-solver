function [M] = map_count_moment(MAP, t, order)
% Computes power moments of counts, at resolution t, of a MAP.
% INPUT
% - MAP: Markovian Arrival Process
% - t: resolution
% - order: orders of the moments to compute
% OUTPUT
% - M: power moments of counts

n = size(MAP{1},1);
theta = map_prob(MAP);

if map_issym(MAP)
    e = sym(ones(n,1));
    M = sym(zeros(length(order),1));
else
    e = ones(n,1);
    M = zeros(length(order),1);
end

if map_issym(MAP) || max(order) > 4
    % symbolic derivative
    z = sym('z');
    MZ = theta*expm(MAP{1}*t+MAP{2}*exp(z)*t)*e;
    for i = 1:length(order)
        M(i) = subs(diff(MZ,z,order(i)),z,0);
    end
else
    % Analytic derivative of M(z) = theta*expm(D0*t + exp(z)*D1*t)*e at z=0,
    % with no step size and no Symbolic Toolbox. Write the exponent as a
    % polynomial in z truncated at z^K,
    %
    %   B(z) = (D0+D1)*t + sum_{j>=1} (D1*t/j!) * z^j,
    %
    % and carry it in the block upper-triangular Toeplitz matrix G whose
    % (p,q) block is B_{q-p}. Such matrices represent polynomials modulo
    % z^(K+1) faithfully, and the representation is a ring homomorphism, so
    % expm(G) represents expm(B(z)): its (1,k+1) block is the z^k Taylor
    % coefficient, that is the k-th derivative over k!. One exponential of
    % size (K+1)*n therefore yields every requested order at once, exact to
    % the tolerance of expm rather than to a differencing step.
    %
    % This replaced the vendored derivest (Romberg extrapolation, orders 1-4
    % only) on 2026-08-19. It is the same method as the C++ twin in
    % cpp/include/line/api/mam/map_count_moment.h; on a Poisson MAP the two
    % now agree with the exact Touchard value at every order, where derivest
    % was 1.7e-5 off at order 4. See _kb/03-api-layer.md.
    K = max(order);
    G = zeros((K+1)*n);
    for p = 0:K
        for q = p:K
            j = q - p;
            if j == 0
                Bj = (MAP{1} + MAP{2}) * t;
            else
                Bj = MAP{2} * (t / factorial(j));
            end
            G(p*n+(1:n), q*n+(1:n)) = Bj;
        end
    end
    EG = expm(G);
    for i = 1:length(order)
        k = order(i);
        M(i) = factorial(k) * theta * EG(1:n, k*n+(1:n)) * e;
    end
end

end
