% https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-10-evaluating-a-ctmc-symbolically
model = gallery_cqn(1,true);

[infGen, eventFilt, syncInfo, stateSpace, nodeStateSpace] = SolverCTMC(model).getSymbolicGenerator;
infGen
stateSpace

syncInfo{1}.active{1}.print, syncInfo{1}.passive{1}.print % x1
syncInfo{2}.active{1}.print, syncInfo{2}.passive{1}.print % x2

steadyStateVector = CTMC(infGen).solve()'

