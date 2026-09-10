function pc = mtrace_pc(~,C)
% Computes the probabilities of arrival for each class.
% T: the inter-arrival times (ignored, for orthogonality with other APIs)
% C: the class labels

labels = unique(C);
m = length(labels);

pc = zeros(m,1);

for i = 1:m
   pc(i) = sum(C==labels(i))/length(C); 
end