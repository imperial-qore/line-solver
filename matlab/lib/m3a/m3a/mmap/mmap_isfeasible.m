function [TF] = mmap_isfeasible(MMAP,TOL)
% Checks whether a MMAP is feasible up to the given tolerance.
% If the tolerance is not specified, the default tolerance mapqntbx_feastol
% is used.

if (nargin == 1)
   TOL = 10^(-mapqntbx_feastol); 
end

for i = 1:length(MMAP)
   if max(max(abs(imag(MMAP{1})))) > TOL
       TF = 0;
       return;
   end
end

% rows of D0 + D1 sum to zero
% diagonal elements of D0 are < 0
% non-diagonal elements of D0 are >= 0
% elements of D1 are >= 0
TF = map_isfeasible(MMAP);
if TF == 0
    return;
end

C = length(MMAP)-2;

% elements of D1c are >= 0
for c = 1:C
   smallest = min(min(MMAP{2+c}));
   if (smallest < -TOL)
       TF = 0;
       return;
   end
end

% D1 = D11 + D12 + ... + D1C
S = MMAP{2};
for c = 1:C
    S = S - MMAP{2+c};
end
if max(max(abs(S))) > TOL
   TF = 0;
   return;
end

end