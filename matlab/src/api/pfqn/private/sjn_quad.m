%{
%{
 % @file sjn_quad.m
 % @brief Quadrature and closed-form integrals for the SJN response time equations.
%}
%}

function varargout = sjn_quad(op, varargin)
%{
%{
 % @brief Dispatcher for the elementary integrals the SJN conditional waiting
 %        time recursion needs. Grouping them in one file keeps the exact
 %        (pfqn_mvasjn) and the approximate (pfqn_amvasjn) recursions on a
 %        single implementation of the quadrature, which must agree for the
 %        two to be comparable. Supported operations:
 %
 %        'pdf'     f(x) of an Erlang mixture on a grid.
 %        'theta'   int_0^x t f(t) dt, closed form.
 %        'ccdf'    int_x^inf f(t) dt, closed form.
 %        'tailmom' int_Lx^inf t^order exp(-c (t-Lx)) f(t) dt, closed form,
 %                  evaluated in logarithms so that exp(c Lx) cannot overflow.
 %        'simpson' composite Simpson over an even number of subdivisions.
 %        'cumsimpson' cumulative Simpson: full panels at the odd nodes and a
 %                  half panel at the even ones, so the primitive is available
 %                  at every grid node, which quadrature at arbitrary abscissae
 %                  could not provide.
 % @fn sjn_quad(op, ...)
 % @param op Operation name, see above.
 % @param varargin Arguments of the selected operation, in the order that
 %          operation's branch reads them.
 % @return varargout Result of the requested operation.
%}
%}
switch op
    case 'pdf'
        varargout{1} = local_pdf(varargin{1}, varargin{2});
    case 'theta'
        varargout{1} = local_theta(varargin{1}, varargin{2});
    case 'ccdf'
        varargout{1} = local_ccdf(varargin{1}, varargin{2});
    case 'tailmom'
        varargout{1} = local_tailmom(varargin{1}, varargin{2}, varargin{3}, varargin{4});
    case 'simpson'
        varargout{1} = local_simpson(varargin{1}, varargin{2});
    case 'cumsimpson'
        varargout{1} = local_cumsimpson(varargin{1}, varargin{2});
    otherwise
        line_error(mfilename,sprintf('unknown operation ''%s''',op));
end
end

function y = local_pdf(f, x)
y = zeros(size(x));
for j = 1:length(f.w)
    k = f.k(j); mu = f.mu(j);
    y = y + f.w(j) * exp(k*log(mu) + (k-1)*log(max(x,realmin)) - mu*x - gammaln(k));
end
end

function y = local_theta(f, x)
y = zeros(size(x));
for j = 1:length(f.w)
    k = f.k(j); mu = f.mu(j);
    y = y + f.w(j) * (k/mu) * gammainc(mu*x, k+1);
end
end

function y = local_ccdf(f, x)
y = 0;
for j = 1:length(f.w)
    y = y + f.w(j) * gammainc(f.mu(j)*x, f.k(j), 'upper');
end
end

function y = local_tailmom(f, Lx, c, order)
y = 0;
for j = 1:length(f.w)
    k = f.k(j); mu = f.mu(j);
    rate = mu + c;
    g = gammainc(rate*Lx, k+order, 'upper');
    if g <= 0
        continue
    end
    lg = c*Lx + k*log(mu/rate) + log(g);
    if order == 1
        lg = lg + log(k/rate);
    end
    y = y + f.w(j) * exp(lg);
end
end

function I = local_simpson(y, dx)
n = numel(y);
I = dx/3 * (y(1) + y(n) + 4*sum(y(2:2:n-1)) + 2*sum(y(3:2:n-2)));
end

function I = local_cumsimpson(y, dx)
n = numel(y);
I = zeros(1,n);
for i = 3:2:n
    I(i) = I(i-2) + dx/3 * (y(i-2) + 4*y(i-1) + y(i));
end
for i = 2:2:n
    if i+1 <= n
        I(i) = I(i-1) + dx/12 * (5*y(i-1) + 8*y(i) - y(i+1));
    else
        I(i) = I(i-1) + dx/12 * (-y(i-2) + 8*y(i-1) + 5*y(i));
    end
end
end
