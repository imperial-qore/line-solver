%  r = CheckMERepresentation(alpha, A, prec)
%  
%  Checks if the given vector and matrix define a valid matrix-
%  exponential representation.
%  
%  Parameters
%  ----------
%  alpha : matrix, shape (1,M)
%      Initial vector of the matrix-exponential distribution 
%      to check
%  A : matrix, shape (M,M)
%      Matrix parameter of the matrix-exponential distribution
%      to check
%  prec : double, optional
%      Numerical precision. The default value is 1e-14.
%  
%  Returns
%  -------
%  r : bool
%      True, if the matrix is a square matrix, the vector and 
%      the matrix have the same size, the dominant eigenvalue
%      is negative and real
%  
%  Notes
%  -----
%  This procedure does not check the positivity of the density!
%  Call 'CheckMEPositiveDensity' if it is needed, but keep in
%  mind that it can be time-consuming, while this procedure
%  is fast.

function r = CheckMERepresentation (alpha, A, prec)

    global BuToolsVerbose;
    global BuToolsCheckPrecision;
    if isempty(BuToolsCheckPrecision)
        BuToolsCheckPrecision = 1e-14;
    end
    
    if ~exist('prec','var')
        prec = BuToolsCheckPrecision;
    end

    if size(A,1)~=size(A,2)
        if BuToolsVerbose
            fprintf ('CheckMERepresentation: The matrix is not a square matrix!\n');
        end
        r = false;
        return;
    end

    if length(alpha)~=size(A,1)
        if BuToolsVerbose
            fprintf ('CheckMERepresentation: The vector and the matrix have different sizes!\n');
        end
        r = false;
        return;
    end

    if sum(alpha)<-prec*length(alpha) || sum(alpha)>1+prec*length(alpha)
        if BuToolsVerbose
            fprintf ('CheckMERepresentation: The sum of the vector elements is less than zero or greater than one (precision: %g)!\n',prec);
        end
        r = false;
        return;
    end

    if max(real(eig(A)))>=prec
        if BuToolsVerbose
            fprintf ('CheckMERepresentation: There is an eigenvalue of the matrix with non-negative real part (at precision %g)!\n',prec);
        end
        r = false;
        return;
    end
    
    ev = eig(A);
    [~,ix] = sort(abs(real(ev)));
    maxev = ev(ix(1));

    % The dominant real part need not be attained by a single eigenvalue, and
    % the sort above breaks such a tie arbitrarily. A concentrated matrix
    % exponential of order 2n+1 puts its whole spectrum on the line
    % Re = -mu1, so the arbitrary pick returns a complex eigenvalue and the
    % test below rejects a valid ME even though the real eigenvalue -mu1 is
    % equally dominant. Among the eigenvalues attaining the dominant real
    % part, prefer a real one: "the dominant eigenvalue is real" is a
    % statement about the spectrum, not about which tied eigenvalue the sort
    % happened to return. This only widens the accepted set, and only in the
    % tied case, which BuTools itself reports as a warning below rather than
    % as a rejection. Matches native Python butools.ph.check and
    % jline.lib.butools.ph.CheckMERepresentation.
    domre = abs(real(maxev));
    domtol = max(prec, 1e-8*domre);
    tied = ev(abs(abs(real(ev)) - domre) <= domtol);
    realtied = tied(abs(imag(tied)) <= domtol);
    if ~isempty(realtied)
        maxev = real(realtied(1));
    end

    % Value-based rather than isreal(maxev): eig returns a complex array as
    % soon as one eigenvalue is complex, so isreal() is false for a real-valued
    % entry of that array and the dominant eigenvalue of any ME with an
    % oscillating component would be declared non-real. jline.lib.butools
    % CheckMERepresentation and native Python already test the imaginary part.
    if abs(imag(maxev)) > prec
        if BuToolsVerbose
            fprintf ('CheckMERepresentation: The dominant eigenvalue of the matrix is not real!\n');
        end
        r = false;
        return;
    end

    % Nested rather than a single && expression: BuToolsVerbose is a global that
    % is [] until lineStart or BuToolsInit runs, and `scalar && []` is a hard
    % error in MATLAB. The left operand here is true exactly when the spectrum
    % has a repeated dominant modulus, which is the common case for an Erlang
    % (ME.fromErlang(2,2) has eigenvalues -2,-2), so the short circuit that hides
    % the problem for a generic ME does not fire and validity checking dies on a
    % verbosity flag. `if []` is simply false, so the nested form is safe
    % whether or not the global has been initialised.
    if sum(abs(ev(1:end))==abs(maxev)) > 1
        if BuToolsVerbose
            fprintf ('CheckMERepresentation warning: There are more than one eigenvalue with the same absolute value as the largest eigenvalue!\n');
        end
    end

    r = true;
end
