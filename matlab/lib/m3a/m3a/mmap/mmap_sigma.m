function sigma = mmap_sigma(mmap)
% Computes one-step class transition probabilities, i.e.
% p_{i,j} = P(C_k = j | C_{k-1} = i)
% INPUT
% - mmap: the MMAP
% OUTPUT
% - sigma: the class-transition probabilities

C = length(mmap)-2;

if map_issym(mmap)
    sigma = sym(zeros(C,C));
else
    sigma = zeros(C,C);
end

alpha = map_pie(mmap);

for i = 1:C
    % alpha * P_i
    start = alpha * (-mmap{1}\mmap{2+i});
    for j = 1:C
        % alpha * P_i * P_j * 1
        sigma(i,j) = sum(start * ((-mmap{1}) \ mmap{2+j}));
    end
end

if map_issym(mmap)
   sigma = simplify(sigma); 
end