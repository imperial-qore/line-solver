% Test gallery_hypm1 with MVA
model = gallery_hypm1();
solver = MVA(model);
avgTable = solver.getAvgTable();
fprintf('Model: %s\n', model.getName());
disp(avgTable);
