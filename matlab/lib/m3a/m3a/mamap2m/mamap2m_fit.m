function MMAP = mamap2m_fit(M1, M2, M3, GAMMA, P, F, B, S, fbsWeights)
% Fits a MAPH(2,m) or A MAMAP(2,m) that matches the provided
% characteristics. Three characteristics of the marked process are matched
% either exactly or approximately. Among these, the class probabilities are
% always matched exactly. The remaining two characteristics, by default,
% are the forward and backward moments, unless the underlying AMAP(2) is
% degenerate. Different pairs of characteristics can be chosen by
% specifying higher weights in the optional parameter fbsWeights.
%
% Input
% - M1,M2,M3: moments of the inter-arrival times
% - GAMMA: auto-correlation decay rate of the inter-arrival times
% - P: class probabilities
% - F: first-order forward moments
% - B: first-order backward moments
% - fbsWeights: the weight assigned to forward moments, backward moments
%               and class transition probabilities (default: [1, 1, 1])
% Output
% - MMAP: fitted MAMAP[m]

if nargin == 2
    % by default equal weights
    % when weights are equal, F+B is preferred over F+S and B+S
    fbsWeights = [1 1 1];
end

fbWeights = [fbsWeights(1), fbsWeights(2)];
fsWeights = [fbsWeights(1), fbsWeights(3)];
bsWeights = [fbsWeights(2), fbsWeights(3)];

gammatol = 1e-4;
degentol = 1e-8;

if length(P) > 2
    fprintf('Fitting MAMAP(2,m): fitting F+B because m > 2\n');
    MMAP = mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B);
    return;
end

% TODO: if it is a Poisson process, convert to a second-order Poisson
% in order to have more degrees of freedom for the marked process

% TODO: if it is a PH-Renewal process and forward moments are preferred,
% convert to non-canonical form

if abs(GAMMA) < gammatol
    fprintf('Fitting MAMAP(2,m): fitting MAPH because gamma = %e\n', GAMMA);
    MMAP = maph2m_fit(M1, M2, M3, P, B);
    return;
end

% fit underlying AMAP(2)
[~,MAPS] = amap2_fit_gamma(M1, M2, M3, GAMMA);

% fit marked Poisson process
if length(MAPS) == 1 && size(MAPS{1}{1},1) == 1
    fprintf('Fitting MAMAP(2,m): fitting marked Poisson because the underlying process has one state\n');
    MAP = MAPS{1};
    % number of classes
    m = length(P);
    % fit class probabilities
    MMAP = cell(1,2+m);
    MMAP{1} = MAP{1};
    MMAP{2} = MAP{2};
    for c = 1:m
        MMAP{2+c} = MMAP{2} .* P(c);
    end
    return;
end

MMAPS = cell(1,length(MAPS));
ERRORS = zeros(1,length(MAPS));
for j = 1:length(MAPS)
    h1 = -1/MAPS{j}{1}(1,1);
    h2 = -1/MAPS{j}{1}(2,2);
    r1 = h1*MAPS{j}{1}(1,2);
    r2 = h2*MAPS{j}{2}(2,2);
    % detect and handle degenerate cases
    degen = 1;
    if GAMMA > 0
        if r1 < degentol || (1-r2) < degentol
            error('Should not happen');
        elseif abs(h2-h1*r2) < degentol
            MMAPS{j} = mamap22_fit_fs_multiclass(MAPS{j}, P, F, S, [], fsWeights);
        elseif abs(h1 - h2 + h2*r1) < degentol
            MMAPS{j} = mamap22_fit_bs_multiclass(MAPS{j}, P, B, S, [], bsWeights);
        elseif (1-r1) < degentol
            % non-canonical APH(2)
            % only forward can be fitted
            MMAPS{j} = mamap22_fit_fs_multiclass(MAPS{j}, P, F, S, [], fsWeights);
        elseif r2 < degentol
            % canonical APH(2)
            MMAPS{j} = maph2m_fit_multiclass(MAPS{j}, P, B);
        else
            degen = 0;
        end
    else
        if (1-r2) < degentol
            error('Should not happen');
        elseif abs(h1 - h2 - h1*r1 + h1*r1*r2) < degentol
            MMAPS{j} = mamap22_fit_fs_multiclass(MAPS{j}, P, F, S, [], fsWeights);
        elseif abs(h1 - h2 + h2*r1) < degentol
            MMAPS{j} = mamap22_fit_bs_multiclass(MAPS{j}, P, B, S, [], bsWeights);
        elseif r2 < degentol && (1-r1) < degentol
            % degenerate canonical APH(2)
            MMAPS{j} = maph2m_fit_multiclass(MAPS{j}, P, B);
        elseif r2 < degentol
            if fbsWeights(1) >= fbsWeights(2)
                % fit forward or sigma
                MMAPS{j} = mamap22_fit_fs_multiclass(MAPS{j}, P, F, S, [], fsWeights);
            else
                % fit backward or sigma
                MMAPS{j} = mamap22_fit_bs_multiclass(MAPS{j}, P, B, S, [], bsWeights);
            end
        else
            degen = 0;
        end
    end
    % handle non-degenerate cases according to user preference
    if ~degen
        if fbsWeights(1) >= fbsWeights(3) && fbsWeights(2) >= fbsWeights(3)
            % prefer forward and backward
            MMAPS{j} = mamap2m_fit_fb_multiclass(MAPS{j}, P, F, B, [], fbWeights);
        elseif fbsWeights(1) >= fbsWeights(2)
            % prefer forward and sigma
            MMAPS{j} = mamap22_fit_fs_multiclass(MAPS{j}, P, F, S, [], fsWeights);
        else
            % prefer backward and sigma
            MMAPS{j} = mamap22_fit_bs_multiclass(MAPS{j}, P, B, S, [], bsWeights);
        end
    end
    % compute fitting error
    fF = mmap_forward_moment(MMAPS{j}, 1);
    fB = mmap_backward_moment(MMAPS{j}, 1);
    fS = mmap_sigma(MMAPS{j});
    ERRORS(j) = fbsWeights(1) * (F(1)/fF(1)-1)^2 + ...
                fbsWeights(2) * (B(1)/fB(1)-1)^2 + ...
                fbsWeights(3) * (S(1,1)/fS(1,1)-1)^2;
end

% pick the best fit
[~,BEST] = min(ERRORS);
MMAP = MMAPS{BEST};

end