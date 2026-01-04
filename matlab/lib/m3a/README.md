M3A: Marked MAP Matching Algorithms
===
<p><b>Codebase</b>
<br>This fork of the main M3A repository offers a customized version of M3A in use within the <a href="https://github.com/imperial-qore/line-solver">LINE solver</>. In particular, the present version has several additional library functions and  does not print out messages on screen.

<p><b>Introduction</b>
<br>Marked Markovian Arrival Processes (MMAPs) are a class of stochastic processes
used to model multi-class correlated arrivals at a queuing system.
M3A is a set of Matlab functions designed for computing the statistical descriptors
of MMAPs and fitting marked traces with MMAPs. The documentation of M3A is available
<a href="https://github.com/Imperial-AESOP/M3A/blob/master/m3a.pdf">here</a>.

<p><b>Installation Requirements</b>
<br>- A recent version of MATLAB, e.g., R2012 or above.
<br>- <a href="https://github.com/kpctoolboxteam/kpc-toolbox/">KPC-Toolbox</a> version 0.3.2 or above.

<p><b>MATLAB Example</b>
<br>% S(i): inter-arrival time between i-th and (i-1)-th arrivals 
<br>% C(i): class of i-th arrival
<br>load example.mat
<br>T = m3afit_init(S,C)
<br>MMAP = m3afit_auto(T,'NumStates',2)

<p><b>License</b>
<br>M3A is released under BSD-3 license.

<p><b>References</b>
<br>
Giuliano Casale, Andrea Sansottera, Paolo Cremonesi,
*Compact Markov-Modulated Models for Multiclass Trace Fitting*,
European J. of Operational Research, INFORMS, 2016.
<a href="http://www.sciencedirect.com/science/article/pii/S0377221716304258">Paper</a>

<p><b>Funding</b>
<br>
The development of this tool has been supported in part by the Horizon 2020 project <a href="http://www.dice-h2020.eu">DICE</a>.
