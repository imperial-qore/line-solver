function X = dmap_sample(DMAP, n)
% Generate n samples from a discrete MAP
% Returns inter-arrival times (integer-valued)
    D0 = DMAP{1};
    D1 = DMAP{2};
    N = size(D0, 1);
    al = dmap_pie(DMAP);
    phase = randsample(N, 1, true, al);
    X = zeros(n, 1);
    for i = 1:n
        t = 0;
        while true
            t = t + 1;
            probs = [D0(phase, :), D1(phase, :)];
            next = randsample(2*N, 1, true, probs);
            if next > N
                phase = next - N;
                X(i) = t;
                break;
            else
                phase = next;
            end
        end
    end
end
