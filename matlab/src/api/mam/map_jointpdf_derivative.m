function gamma = map_jointpdf_derivative(MAP, iset)
% partial derivative at 0 of a MAP's joint PDF
% A. Horvath et al. A Joint Moments Based Analysis of Networks of
% MAP/MAP/1 Queues

n = length(MAP{1});

gamma = map_pie(MAP);
for j=iset(:)'
    gamma = gamma * MAP{1}^j * MAP{2};
end
gamma = gamma * ones(n,1);
end