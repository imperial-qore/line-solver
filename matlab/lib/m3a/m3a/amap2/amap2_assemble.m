function AMAP = amap2_assemble(l1,l2,p1,p2,form)
% Returns an AMAP(2) with the given parameters.

if form == 1
   AMAP = {[-1/l1  p1/l1; 0  -1/l2], ...
           [(1-p1)/l1  0; (1-p2)/l2  p2/l2]}; 
elseif form == 2
   AMAP = {[-1/l1  p1/l1; 0  -1/l2], ...
           [0  (1-p1)/l1; (1-p2)/l2  p2/l2]};
else
    error('Invalid form: should be either 1 (gamma > 0) or 2 (gamma < 0)');
end