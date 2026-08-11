function varargout = qbd_fundmat(B, L, F, matrices, precision, maxNumIt)
%QBD_FUNDMAT  Cyclic-reduction fundamental matrices G (and optionally R) of a QBD.
%   Operates directly on the raw level blocks (B,L,F) of a homogeneous QBD and
%   supports complex-valued L (as required when L is shifted by -s*I in a
%   Laplace-domain transient analysis). This differs from QBD_RG, which builds
%   the blocks from MAP representations and is only exercised at real argument.
%
%   VARARGOUT = QBD_FUNDMAT(B, L, F, MATRICES) returns the matrices named by
%   the character codes in MATRICES (default 'G'), in that order. Supported
%   codes: 'G', 'R', 'GR', 'RG'.
%
%   The G matrix is obtained by cyclic reduction (Bini-Meini logarithmic
%   reduction); R is recovered from G as R = Fm*(I-(Lm+Fm*G))^-1.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
if nargin < 4 || isempty(matrices);  matrices  = 'G';    end
if nargin < 5 || isempty(precision); precision = 1e-14;  end
if nargin < 6 || isempty(maxNumIt);  maxNumIt  = 50;     end

m = size(L, 1);
II = eye(m);
lamb = max(-real(diag(L)));
Bm = B / lamb;
Lm = L / lamb + II;
Fm = F / lamb;

BF = (II - Lm) \ II;
BB = BF * Fm;
BF = BF * Bm;
G  = BF;
PI = BB;
check = 1;
numit = 0;
while check > precision && numit < maxNumIt
    Lstar = BF * BB + BB * BF;
    Bstar = BB * BB;
    Fstar = BF * BF;
    BB = (II - Lstar) \ II;
    BF = BB * Fstar;
    BB = BB * Bstar;
    G  = G  + PI * BF;
    PI = PI * BB;
    check = min(norm(BB, inf), norm(BF, inf));
    numit = numit + 1;
end

R = [];
outs = cell(1, length(matrices));
for i = 1:length(matrices)
    c = matrices(i);
    if c == 'G'
        outs{i} = G;
    elseif c == 'R'
        if isempty(R)
            R = Fm * ((II - (Lm + Fm * G)) \ II);
        end
        outs{i} = R;
    else
        line_error(mfilename, sprintf('unknown matrix code ''%c'' in ''%s''', c, matrices));
    end
end
varargout = outs;
end
