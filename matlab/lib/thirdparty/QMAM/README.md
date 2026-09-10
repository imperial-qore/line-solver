# Q-MAM #

Q-MAM is a set of MATLAB scripts for the analysis of queues using matrix-analytic methods. 

It includes, among others, the following queueing models:

* PH/PH/1
* MAP/MAP/1
* MAP/M/c 
* MAP/D/c
* RAP/RAP/1
* MMAP[K]/PH[K]/1
* MMAP[K]/SM[K]/1
* SM[K]/PH[K]/1

State-of-the-art solution techniques are used to solve these models efficiently. 


### Quick start ###

* The easiest way to start is to go to the [ Downloads ](https://bitbucket.org/qmam/qmam/downloads) section and download the latest version of QMAM (currently 1.1). 
* The .zip contains a set of matlab scripts that you can add to your path and use directly. 
* Download and add to your path the SMC Solver scripts, from the sources indicated below in the Requirements section. 
* Download [ this example ](https://bitbucket.org/qmam/qmam/src/8df8896007438d42f6b3da8359a878be18e5a8fc/examples/QMAM_example_1.m?at=master) and run it in MATLAB. The examples shows how to use the QMAM scripts to obtain the queue-length and response-time distributions for a PH/PH/1 queue. 
* Detlais on the parameters of each script are provided in the respective script. 

### Requirements ###
Q-MAM requires the [ SMC Solver ](http://win.uantwerpen.be/~vanhoudt/)

* Download the SMC Solver QBD files from [ here ](http://win.uantwerpen.be/~vanhoudt/tools/QBDfiles.zip)
* Download the SMC Solver MG1 files from [ here ](http://win.uantwerpen.be/~vanhoudt/tools/MG1files.zip) 

### License ###

Q-MAM is released under the [ Apache License 2.0 ](http://www.apache.org/licenses/LICENSE-2.0)
### Citing Q-MAM ###

When using this tool, please refer to the paper 

* Q-MAM: A Tool for Solving Infinite Queues using Matrix-Analytic Methods, by J. F. Pérez, J. Van Velthoven and B. Van Houdt. 
Proceedings of Valuetools 2008. [ ACM DL ](http://dl.acm.org/citation.cfm?id=1536977&dl=ACM&coll=DL&CFID=482820390&CFTOKEN=54023148)


Additional info on the tool is found in the paper.

### Support ###
* To report bugs and request features please use the Issue tracker on the left-hand panel.
* For questions about using Q-MAM please contact the repo admin. 

### Contributors ###
* Juan F. Pérez
* Benny Van Houdt
* Jeroen Van Velthoven