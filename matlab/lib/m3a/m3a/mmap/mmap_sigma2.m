function sigma = mmap_sigma2(mmap)
% Computes two-step class transition probabilities, i.e.
% p_{i,j,h} = P(C_k = h | C_{k-1} = j | C_{k-2} = i)
% INPUT
% - mmap: the MMAP
% OUTPUT
% - sigma: the 3D matrix of class-transition probabilities

C = length(mmap)-2;

if map_issym(mmap)
    sigma = sym(zeros(C,C,C));
else
    sigma = zeros(C,C,C);
end

alpha = map_pie(mmap);

for i = 1:C
    starti = alpha * (-mmap{1}\mmap{2+i});
    for j = 1:C
        startj = starti * ((-mmap{1}) \ mmap{2+j});
        for h = 1:C
            sigma(i,j,h) = sum(startj * ((-mmap{1}) \ mmap{2+h}));
        end
    end
end

if map_issym(mmap)
   sigma = simplify(sigma); 
end