function nu = map_ccdf_derivative(MAP, i)
% derivative at 0 of a MAP CCDF
% A. Horvath et al. A Joint Moments Based Analysis of Networks of
% MAP/MAP/1 Queues

n = length(MAP{1});

pie = map_pie(MAP);
nu = pie * MAP{1}^i *ones(n,1);
end