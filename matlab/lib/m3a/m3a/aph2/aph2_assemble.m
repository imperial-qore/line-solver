function APH = aph2_assemble(l1, l2, p1)
% Given the parameters, return the APH(2).

APH = {[-1/l1 1/l1*p1; 0 -1/l2], ...
       [1/l1*(1-p1) 0; 1/l2 0] };

end