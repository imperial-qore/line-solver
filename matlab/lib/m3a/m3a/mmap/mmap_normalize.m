function MMAP = mmap_normalize(MMAP)
% Fixes MMAP feasibility by setting negative values to zero and forcing
% the other conditions.

if isempty(MMAP)    
    return
end
K = size(MMAP{1},1);
C = length(MMAP)-2;

for i = 1:K
   for j = 1:K
       if i ~= j
           MMAP{1}(i,j) = max(MMAP{1}(i,j), 0);
       end
   end
end

MMAP{2} = 0 *MMAP{1};
for c = 1:C
   MMAP{2+c}(MMAP{2+c} < 0) = 0; 
   if isnan(MMAP{2+c})
       MMAP{2+c} = zeros(size(MMAP{2+c}));
   end
   MMAP{2} = MMAP{2} + MMAP{2+c};
end

for k = 1:K
   MMAP{1}(k,k) =  0;
   MMAP{1}(k,k) = -sum(MMAP{1}(k,:)) - sum(MMAP{2}(k,:));
end
end