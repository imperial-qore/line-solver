function al = dmap_pie(DMAP)
% Stationary vector at arrival epochs for a discrete MAP
% DMAP is a cell array {D0, D1}
    D0 = DMAP{1};
    D1 = DMAP{2};
    N = size(D0, 1);
    P = inv(eye(N) - D0) * D1;
    al = dtmc_solve(P);
end
