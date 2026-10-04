# rmftool - A library to Compute (Refined) Mean Field Approximations


This python library implements algorithms to simulate and study population processes, to compute their mean-field approximations and refined mean-field approximations for the transient regime and steady-state.


The tool accepts three model types:
* homogeneous population processes (HomPP)
* density dependent population processes (DDPPs)
* heterogeneous population models (HetPP)

In particular, it provides a numerical algorithm to compute the constant of the "refined mean field approximation" provided in the paper [A Refined Mean Field Approximation](https://hal.inria.fr/hal-01622054/document) by N. Gast and B. Van Houdt, accepted at SIGMETRICS 2018. And a framework to compute heterogeneous mean field approximations [Mean Field and Refined Mean Field Approximations for Heterogeneous Systems: It Works!](https://arxiv.org/pdf/2111.01594.pdf) by N. Gast and S. Allmeier, accepted at SIGMETRICS 2022.

## Installation

### Using Pip

```
pip install rmftool
```

### Using Github

Clone the repository. For example by using the termin git client:
```
git clone https://github.com/ngast/rmf_tool
```
then manually add the package by using pip in the corresponding folder 
```
pip install .
```

## Documentation / Usage

An extensive description of the supported classes can be found in the publication [rmf tool – A library to Compute (Refined) Mean Field Approximation(s)](https://www.performance2021.deib.polimi.it/wp-content/uploads/2021/11/TOSME21_paper_2.pdf) by N. Gast and S. Allmeier. The example given in the publication can be found in the repository [https://gitlab.inria.fr/gast/toolpaper_rmf](https://gitlab.inria.fr/gast/toolpaper_rmf)

### Homogeneous Population Process

Homogeneous population model consists of $n$ interacting objects that each evolve in a finite state space $\{1...d\}$. All objects have similar transition rates which are a combination of unilateral and pairwise interactions, i.e. objects can change their state with or without interacting with one other member of the population. We assume that $X = (X_1...X_d)$ is a continuous time Markov chain whose transitions are such that for all state $s,s′, \tilde{s}, \tilde{s}′$:

- (Unilateral) An object in state $s$ moves to state $s'$ at rate $a_{s,s'}$.
- (Pairwise) A pair of objects in state $(s,\tilde{s})$ moves to state $(s',\tilde{s}')$ at rate $b_{s,\tilde{s},s',\tilde{s}'}/n$.

Homogeneous population models is a subclass of density dependent population processes
(DDPPs) defined in the following. 

### Density Dependent Population Process


A density dependent process is a Markov chain that evolves in a sub-domain of $R^d$ where the transitions are given by a set $L$ and a list of rate-functions $\beta_\ell$. To a system with size $N$ a Markov chain is associated whose transitions are (for all $\ell\in L$):

* $x \mapsto x + \frac 1N \ell$ at rate $\beta_\ell(x)$

The class DDPP can be used to defined a process, via the functions ```add_transition(l, beta)```.


### Heterogeneous Population Process

As before, the heterogeneous population model consists of $n$ interacting objects which each evolve in a finite state space $\{1...d\}$. Each object has a specific transition rate which is a combination of unilateral and pairwise interactions. In contrast to the HomPP, transition rates are object dependent:

- The object $k$ moves from state $s$ to state $s'$ at rate $a_{k,s,s'}$.
- The pair $(k,k')$ moves from states $(s,\tilde{s})$ to states $(s',\tilde{s}')$ at rate $b_{k,k',s,\tilde{s},s',\tilde{s}'}/n$.

------

### (Some) Available Methods

In order to obtain sample trajectories, mean field and refined mean field approximations the following methods can be used.

| Description                         | Function DDPP \& HomPP                               | Function HetPP                 |
| ----------------------------------- | -------------------------------------------- | --------------------------------------- |
| sample trajectoy                    | ```simulate(N, time)```                      | ``` simulate(time)```                   |
| mean-field (or fluid) approximation | ```ode(time)```                              | ```ode(time)```                         |
| refined mean-field  (transient)     | ```meanFieldExpansionTransient(time)```      | ```meanFieldExpansionTransient(time)``` |
| refined mean-field (steady-state)   | ```meanFieldExpansionSteadyState()```        | ```meanFieldExpansionSteadyState()```   |
| steady-state using simulation       | ```steady_state_simulation(N, time=20000)``` | -                                       |
| sample Mean trajectory              | -                                            | ```sampleMean(time, samples=1000)```    |


The functions are documented and their documentation is accessible by the "help" command.

Apart from that, the **documentation** is mostly **contained in the [examples](#examples)** below (from basic to more advanced). 



### Computation time 

#### DDPP

The computation of the function ```meanFieldExpansionSteadyState``` grows as $d^3$, where d is the dimension of the model. When using the option ```symbolic_differentiation=True``` (default), it takes about 10 seconds for a 50-dimensional system, about 60 seconds for a 100-dimensional system. Note that for large models, most of the time is taken by the computation of the symbolic derivatives. This time can be reduced when using ```symbolic_differentiation=False``` to approximatively 4 sec for 50-dimensional systems and about 15 seconds for 100-dimensional systems (note: due to limitation of the lambdify function, the symbolic_differentiation=False cannot be used for systems with more than 255 dimensions). 

-----

The simulation of the underlying Markov chain is not optimized and therefore might be slow for large models. 

## Examples

### Introductory example 

The following code illustrates how to define a 2-dimensional model using the [DDPP](#density-dependent-population-process) class. It plots the mean-field approximation versus one sample trajectory. It then computes the refined mean-field approximation (in steady-state)

```
import rmftool as rmf

ddpp = rmf.DDPP()

ddpp.add_transition([-1,1], lambda x : x[0])
ddpp.add_transition([1,-1], lambda x : x[1]+x[0]*x[1])

ddpp.set_initial_state([.5,.5])

ddpp.plot_ODE_vs_simulation(N=100)
```

### SIS model 

One of the simplest examples of population process is the epidemic model called the SIS model. In an SIS model, each object can be in one of the two states $S$ (susceptible) or $I$ (infected). Susceptible objects can become infected from an external source (unilateral transition) or when meeting an infected individual (pairwise transition). An infected individual can recover and become susceptible again (unilateral transition). We assume that an individual becomes infected at rate $\alpha$ by an external source, and recovers at rate $\beta$. Moreover, assume that the rate at which two individuals meet is $\gamma/n$ and that when a susceptible meets an infected individual, the susceptible gets infected. 

In the following, we illustrate the setup using the different available classes.



With our tool, we define a class called HomPP for which we specify the transition rates and an initial state. For the SIS model above, with $\alpha, \beta, \gamma = 1$, we write:

#### Homogeneous Population Process (HomPP)

```
import rmftool as rmf
import numpy as np

model = rmf.HomPP()

d, S, I = 2, 0, 1

A, B = np.zeros((d, d)), np.zeros((d, d, d, d)) 

A[S, I] = 1                 # unil. transitions from S to I (alpha)
A[I, S] = 1                 # unil. transitions from I to S (beta)
B[S, I, I, I] = 1           # pairwise transition from S to I (gamma)

model.add_rate_tensors(A, B)

model.set_initial_state([1,0])
```


#### Density Dependent Population Process (DDPP)

Using the class DDPP we can define density dependent population processes directly from their mathematical definition. For the above SIS example, we write: 

```
import rmftool as rmf

model = rmf.DDPP()

alpha, beta, gamma = 1, 1, 1

model.add_transition([-1,1], lambda x: alpha*x[0])
model.add_transition([1,-1], lambda x: beta*x[1])
model.add_transition([-1,1], lambda x: gamma*x[0]*x[1])

model.set_initial_state([1,0])
```

#### Heterogeneous Population Processs (HetPP)

To set up a heterogeneous version of the previous SIS model we use the HetPP class of the toolbox. In contrast to the HomPP and DDPP class, the model **can not** be defined independent of the system size, i.e., $n$ and $d$ have to be defined to initialize the model. For instance, to set up a SIS model where objects are almost identical but some recover slower than others, we can use the code:

```
import rmftool as rmf
import numpy as np

model = rmf.HetPP()

N, d = 20, 2
S, I = 0, 1

A, B = np.zeros((N, d, d)), np.zeros((N, N, d, d, d, d))

A[:, S, I] = np.ones((N))
A[:, I, S] = np.random.rand(N)   # Hetero. recovery rates
B[:, :, S, I, I, I] = (1/N) * np.ones((N, N))

model.add_rate_tensors(A, B)

model.set_initial_state(np.ones((N,d))*np.array([1,0]))
```

------

Following the setup the methods described in the subsection ['available functions'](#some-available-methods) can be used to analyze the model. We point to 

* [Simple SIR model](examples/BasicExample_SIR.ipynb)

for a more detailed discussion on how to use the toolbox. 

### Advanced examples

* [2-choice model](examples/Example_2choice.ipynb) 
* [Non stable SIR](examples/Example_nonStableSIR.ipynb)
* [Bike sharing system](examples/Example_bikeSharingSystem.ipynb)

More examples can be found in the [example directory](examples) and in the rmf_tool paper repository - [https://gitlab.inria.fr/gast/toolpaper_rmf/-/tree/master/code](https://gitlab.inria.fr/gast/toolpaper_rmf/-/tree/master/code).

## Dependencies

This library depends on the following python library:

* numpy
* random
* scipy
* sympy 
* matplotlib.pyplot
* jax

## History
* v0.3:
  - (New feature) Added two new population classes (HomPP, HetPP) with support for simulation, mean field and first order refined mean field approximation 
* v0.2:
  - (New feature) Support for multi-class variables
  - (Optimization) Addition of an option to use a numerical differentiation
* v0.1: original version 

**Copyright** 2017. **Authors** : Nicolas Gast, Emmanuel Rodriguez, Sebastian Allmeier
