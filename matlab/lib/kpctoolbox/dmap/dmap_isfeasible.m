function bool = dmap_isfeasible(DMAP)
% Check if the DMAP representation is feasible
% DMAP is a cell array {D0, D1}
    D0 = DMAP{1};
    D1 = DMAP{2};
    bool = true;
    if size(D0,1) ~= size(D0,2) || size(D1,1) ~= size(D1,2)
        bool = false; return;
    end
    if size(D0,1) ~= size(D1,1)
        bool = false; return;
    end
    if any(D0(:) < -1e-10) || any(D1(:) < -1e-10)
        bool = false; return;
    end
    P = D0 + D1;
    rs = sum(P, 2);
    if any(abs(rs - 1) > 1e-6)
        bool = false; return;
    end
    d1rs = sum(D1, 2);
    if any(d1rs < -1e-10)
        bool = false; return;
    end
end
