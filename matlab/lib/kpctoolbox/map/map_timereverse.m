function MAPr=map_timereverse(MAP)
piq = map_prob(MAP);
D=diag(piq);
MAPr={inv(D)*MAP{1}'*D,inv(D)*MAP{2}'*D};
end
