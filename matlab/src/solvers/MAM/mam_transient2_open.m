function V = mam_transient2_open(B, L, F, Lv, T, n, m, s)
%MAM_TRANSIENT2_OPEN  Laplace-domain transient V(s,n,m) for an open (infinite) QBD.
%   V = MAM_TRANSIENT2_OPEN(B, L, F, Lv, T, n, m, s) returns the Laplace
%   transform (at complex argument s) of the transient transition-probability
%   matrix from level n to level m of a piecewise level-dependent QBD with
%   regime thresholds T (the last regime repeats to infinity).
%
%   Block cell arrays are indexed by regime k = 1..K (K = length(T)):
%     B{k}  backward (level-down) block for regime k
%     L{k}  local (level-internal) block used on repeating levels of regime k
%     F{k}  forward (level-up) block for regime k
%     Lv{k} local block on the boundary level T(k) of regime k
%
%   Ported from the transient-QBD research code (Horvath et al. formulation).
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
    K = length(T);

    Gs  = cell(K, 1);
    Rs  = cell(K, 1);
    Ghs = cell(K, 1);
    Rhs = cell(K, 1);
    for k = 1:K
        if k < K && T(k+1) - T(k) == 1
            Gs{k}  = [];
            Rs{k}  = [];
            Ghs{k} = [];
            Rhs{k} = [];
        else
            Ik = eye(size(L{k}, 1));
            [G, R] = qbd_fundmat(B{k}, L{k} - s*Ik, F{k}, 'GR');
            Gs{k} = G;  Rs{k} = R;
            [G, R] = qbd_fundmat(F{k}, L{k} - s*Ik, B{k}, 'GR');
            Ghs{k} = G; Rhs{k} = R;
        end
    end

    SvHn  = cell(K-1, 1);
    SvH0  = cell(K-1, 1);
    SvHhn = cell(K-1, 1);
    SvHh0 = cell(K-1, 1);
    for k = 1:K-1
        NN = size(Lv{k}, 1);
        if T(k+1) - T(k) > 1
            d  = T(k+1) - T(k);
            Ik = eye(NN);
            num = [mxpow(Ghs{k}, d-1), Gs{k};  Ghs{k}, mxpow(Gs{k}, d-1)];
            den = [Ik,               mxpow(Gs{k}, d);
                   mxpow(Ghs{k}, d), Ik];
            SH = num / den;
            SvHn{k}  = SH(1:NN,      1:NN);
            SvH0{k}  = SH(1:NN,      NN+1:2*NN);
            SvHhn{k} = SH(NN+1:2*NN, 1:NN);
            SvHh0{k} = SH(NN+1:2*NN, NN+1:2*NN);
        else
            NN1 = size(Lv{k+1}, 1);
            SvH0{k}  = zeros(NN1, NN);
            SvHh0{k} = eye(NN);
            SvHhn{k} = zeros(NN, NN1);
            SvHn{k}  = eye(NN1);
        end
    end

    SY = cell(K, 1);
    SY{K} = Gs{K};
    for k = K-1:-1:1
        NNk1 = size(Lv{k+1}, 1);
        SY{k} = SvH0{k} + SvHn{k} * ((s*eye(NNk1) - Lv{k+1} - F{k+1}*SY{k+1} - B{k}*SvHhn{k}) \ (B{k} * SvHh0{k}));
    end

    NN1 = size(Lv{1}, 1);
    SYh = cell(K-1, 1);
    if K >= 1
        SYh{1} = SvHhn{1} + SvHh0{1} * ((s*eye(NN1) - Lv{1} - F{1}*SvH0{1}) \ (F{1} * SvHn{1}));
    end
    for k = 2:K-1
        NNk = size(Lv{k}, 1);
        SYh{k} = SvHhn{k} + SvHh0{k} * ((s*eye(NNk) - Lv{k} - B{k-1}*SYh{k-1} - F{k}*SvH0{k}) \ (F{k} * SvHn{k}));
    end

    SV = cell(K, K);
    for l = 0:K-1
        if l == 0
            SV{l+1, l+1} = (s*eye(NN1) - Lv{1} - F{1}*SY{1}) \ eye(NN1);
        else
            NNl1 = size(Lv{l+1}, 1);
            SV{l+1, l+1} = (s*eye(NNl1) - Lv{l+1} - F{l+1}*SY{l+1} - B{l}*SYh{l}) \ eye(NNl1);
        end
        for k = l+1:K-1
            NNk1 = size(Lv{k+1}, 1);
            SV{k+1, l+1} = (s*eye(NNk1) - Lv{k+1} - F{k+1}*SY{k+1} - B{k}*SvHhn{k}) \ (B{k} * SvHh0{k} * SV{k, l+1});
        end
        for k = l-1:-1:1
            NNk1 = size(Lv{k+1}, 1);
            SV{k+1, l+1} = (s*eye(NNk1) - Lv{k+1} - F{k+1}*SvH0{k+1} - B{k}*SYh{k}) \ (F{k+1} * SvHn{k+1} * SV{k+2, l+1});
        end
        if l > 0
            SV{1, l+1} = (s*eye(NN1) - Lv{1} - F{1}*SvH0{1}) \ (F{1} * SvHn{1} * SV{2, l+1});
        end
    end

    pos = find(T > n, 1);
    if isempty(pos), kn = K; else, kn = pos - 1; end
    pos = find(T > m, 1);
    if isempty(pos), km = K; else, km = pos - 1; end

    NN = size(Lv{kn}, 1);
    II = eye(NN);

    Vu = []; Lu = [];
    thresholdN = (T(kn) == n);
    if thresholdN
        if T(km) == m
            V = SV{kn, km};
            return;
        end
        Vl = SV{kn, km};
        Ll = T(km);
        if km < K
            Vu = SV{kn, km+1};
            Lu = T(km+1);
        end
    else
        if kn < K
            d1 = T(kn+1) - n;
            num = [mxpow(Ghs{kn}, d1-1), Gs{kn};  Ghs{kn}, mxpow(Gs{kn}, d1-1)];
            den = [II, mxpow(Gs{kn}, d1); mxpow(Ghs{kn}, d1), II];
            Tmp = num / den;
            HTnn  = Tmp(1:NN,      1:NN);
            HTn0  = Tmp(1:NN,      NN+1:2*NN);
            HhTnn = Tmp(NN+1:2*NN, 1:NN);
            HhTn0 = Tmp(NN+1:2*NN, NN+1:2*NN);

            d2 = n - T(kn);
            num = [mxpow(Ghs{kn}, d2-1), Gs{kn};  Ghs{kn}, mxpow(Gs{kn}, d2-1)];
            den = [II, mxpow(Gs{kn}, d2); mxpow(Ghs{kn}, d2), II];
            Tmp = num / den;
            HnTn  = Tmp(1:NN,      1:NN);
            HnT0  = Tmp(1:NN,      NN+1:2*NN);
            HhnTn = Tmp(NN+1:2*NN, 1:NN);
            HhnT0 = Tmp(NN+1:2*NN, NN+1:2*NN);

            NNkn1 = size(Lv{kn+1}, 1);
            Yn = HTn0 + HTnn * ((s*eye(NNkn1) - Lv{kn+1} - F{kn+1}*SY{kn+1} - B{kn}*HhTnn) \ (B{kn} * HhTn0));
        else
            d2 = n - T(kn);
            num = [mxpow(Ghs{kn}, d2-1), Gs{kn};  Ghs{kn}, mxpow(Gs{kn}, d2-1)];
            den = [II, mxpow(Gs{kn}, d2); mxpow(Ghs{kn}, d2), II];
            Tmp = num / den;
            HnTn  = Tmp(1:NN,      1:NN);
            HnT0  = Tmp(1:NN,      NN+1:2*NN);
            HhnTn = Tmp(NN+1:2*NN, 1:NN);
            HhnT0 = Tmp(NN+1:2*NN, NN+1:2*NN);
            HTnn = []; HTn0 = []; HhTnn = []; HhTn0 = [];
            Yn = Gs{kn};
        end

        if kn == 1
            Yhn = HhnTn + HhnT0 * ((s*eye(NN1) - Lv{1} - F{1}*HnT0) \ (F{1} * HnTn));
        else
            NNkn = size(Lv{kn}, 1);
            Yhn = HhnTn + HhnT0 * ((s*eye(NNkn) - Lv{kn} - B{kn-1}*SYh{kn-1} - F{kn}*HnT0) \ (F{kn} * HnTn));
        end

        Mkn = s*eye(size(L{kn},1)) - L{kn};
        if T(km) < n
            Vnl = (Mkn - B{kn}*HhnTn - F{kn}*Yn) \ (B{kn} * HhnT0 * SV{kn, km});
        else
            Vnl = (Mkn - F{kn}*HTn0 - B{kn}*Yhn) \ (F{kn} * HTnn * SV{kn+1, km});
        end
        if m == T(km)
            V = Vnl;
            return;
        end

        Vnu = [];
        if km < K
            if T(km+1) < n
                Vnu = (Mkn - B{kn}*HhnTn - F{kn}*Yn) \ (B{kn} * HhnT0 * SV{kn, km+1});
            else
                Vnu = (Mkn - F{kn}*HTn0 - B{kn}*Yhn) \ (F{kn} * HTnn * SV{kn+1, km+1});
            end
        end
        Vnn = (s*II - L{kn} - B{kn}*Yhn - F{kn}*Yn) \ II;
        if n == m
            V = Vnn;
            return;
        end
        if km == K && n < m
            if kn < K
                V = Vnl * mxpow(Rs{km}, m - T(km));
            else
                V = Vnn * mxpow(Rs{km}, m - n);
            end
            return;
        end
        if kn ~= km
            Vu = Vnu; Vl = Vnl; Lu = T(km+1); Ll = T(km);
        elseif n <= m
            Vu = Vnu; Vl = Vnn; Lu = T(km+1); Ll = n;
        else
            Vu = Vnn; Vl = Vnl; Lu = n;       Ll = T(km);
        end
    end

    if km == K && n < m
        if kn < K
            V = Vl * mxpow(Rs{km}, m - T(km));
        else
            V = SV{kn, kn} * mxpow(Rs{km}, m - n);
        end
        return;
    end

    NN = size(Rs{km}, 1);
    II = eye(NN);
    Zden = [II, mxpow(Rs{km}, Lu-Ll); mxpow(Rhs{km}, Lu-Ll), II];
    Znum = [mxpow(Rs{km}, m-Ll); mxpow(Rhs{km}, Lu-m)];
    Z = Zden \ Znum;
    V = Vl * Z(1:NN, :) + Vu * Z(NN+1:2*NN, :);
end
