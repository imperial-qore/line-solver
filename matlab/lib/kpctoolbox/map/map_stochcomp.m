function [MAPnew] = map_stochcomp(MAP, retain_idx)
D0 = MAP{1};
D1 = MAP{2};
Q = D0 + D1;

n = size(Q, 1);
eliminated_idx = setdiff(1:n, retain_idx); % eliminated state indices

Q_RE = Q(retain_idx, eliminated_idx);
Q_EE = Q(eliminated_idx, eliminated_idx);
Q_RR = Q(retain_idx, retain_idx);

QNew = Q_RR + Q_RE * (-Q_EE \ Q(eliminated_idx, retain_idx));

D0new = QNew - D1(retain_idx, retain_idx);

D1new = D1(retain_idx, retain_idx) + ...
        Q_RE * (-Q_EE \ D1(eliminated_idx, retain_idx));

MAPnew = map_normalize({D0new, D1new});
end
