clear; 
evalc('cdf_respt_closed');
exampleName = 'cdf_respt_closed';
TOL=1e-2; % 1-percent tolerance
% ONE GATE COVERS ALL FOUR SOURCES. The fluid rows are a quadrature over a
% curve the integrator laid down, so unlike the theory rows (closed form) and
% the simulated rows (a fixed seed) they move with the hook in
% options.odesolvers.accurateStiffOdeSolver. Measured, the two hooks differ by
% 3.5e-03 on the SCV row and 5.1e-04 on the mean, so 1 percent clears both --
% but only under the midpoint rule below, which is what shrank that spread
% from 1.0e-02.

AvgRespTfromTheory_saved = [
    0.1000
    1.0000];
try
    assert(max(max(abs(AvgRespTfromTheory-AvgRespTfromTheory_saved)))<TOL,sprintf('%s changed on AvgRespTfromTheory.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(AvgRespTfromTheory),mat2str(AvgRespTfromTheory_saved))
end

SqCoeffOfVariationRespTfromTheory_saved = [
   1.000000000000000
   0.333333333333333];
try
    assert(max(max(abs(SqCoeffOfVariationRespTfromTheory-SqCoeffOfVariationRespTfromTheory_saved)))<TOL,sprintf('%s changed on SqCoeffOfVariationRespTfromTheory.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(SqCoeffOfVariationRespTfromTheory),mat2str(SqCoeffOfVariationRespTfromTheory_saved))
end

% THE FLUID ROWS ARE READ WITH A MIDPOINT RIEMANN-STIELTJES SUM since
% 2026-08-29, and re-recorded under it. The example used to charge each
% interval's mass dF at the interval's RIGHT endpoint for the fluid curve as
% well as the simulated one. That is exact for a SIMULATED CDF, which jumps AT
% its samples, and first-order for a FLUID one, which is a continuous law read
% off an integrator grid -- so the fluid rows carried a bias that moved with
% whatever grid the integrator chose. It moved: this row read 0.100603 in the
% 2026-08-12 recording, 0.100487 under @ode15s on 2026-08-28, and the SCV row
% below tripped its tolerance outright when the default hook became
% @lsoda_accurate_stiff. Under the midpoint rule the delay row reads 0.0998967
% and the queue row 0.999905, against the exact 0.1 and 1 of the theory row
% above, and the two integrators now agree to 5.1e-04 instead of 5.1e-03.
AvgRespTfromCDFFluid_saved = [
   0.0998966994046832
   0.999904978544848];
try
    assert(max(max(abs(AvgRespTfromCDFFluid-AvgRespTfromCDFFluid_saved)))<TOL,sprintf('%s changed on AvgRespTfromCDFFluid.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(AvgRespTfromCDFFluid),mat2str(AvgRespTfromCDFFluid_saved))
end

% Re-recorded 2026-08-29 alongside the mean above, and under the same midpoint
% rule. This is the row the quadrature was costing the most: the right-endpoint
% sum read 1.00699 under @ode15s and 1.01715 under @lsoda_accurate_stiff, a
% 1.0e-02 spread that a 1% gate cannot straddle. The midpoint sum reads 1.00032
% and 1.00384 on the same two hooks, against the exact 1.0 of the theory row --
% so BOTH integrators now clear the 1 percent gate, and the row is accurate
% rather than merely reproducible.
SqCoeffOfVariationRespTfromCDFFluid_saved = [
   1.00032205287462
   0.333010606979001];
try
    assert(max(max(abs(SqCoeffOfVariationRespTfromCDFFluid-SqCoeffOfVariationRespTfromCDFFluid_saved)))<TOL,sprintf('%s changed on SqCoeffOfVariationRespTfromCDFFluid.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(SqCoeffOfVariationRespTfromCDFFluid),mat2str(SqCoeffOfVariationRespTfromCDFFluid_saved))
end

AvgRespTfromCDFSim_saved = [
   0.099997995814703
   0.997706224358445];
try
    assert(max(max(abs(AvgRespTfromCDFSim-AvgRespTfromCDFSim_saved)))<TOL,sprintf('%s changed on AvgRespTfromCDFSim.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(AvgRespTfromCDFSim),mat2str(AvgRespTfromCDFSim_saved))
end


SqCoeffOfVariationRespTfromCDFSim_saved = [
   1.011017451784501
   0.336453886527369];
try
    assert(max(max(abs(SqCoeffOfVariationRespTfromCDFSim-SqCoeffOfVariationRespTfromCDFSim_saved)))<TOL,sprintf('%s changed on SqCoeffOfVariationRespTfromCDFSim.',exampleName));
catch me
    fprintf('Assertion failed in %s.\n%s\nNew: %s\nOld: %s\n',exampleName,me.message,mat2str(SqCoeffOfVariationRespTfromCDFSim),mat2str(SqCoeffOfVariationRespTfromCDFSim_saved))
end