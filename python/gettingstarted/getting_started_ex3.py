#!/usr/bin/env python
# coding: utf-8

# https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-3-machine-interference-problem

# In[1]:


from line_solver import *


# In[2]:


model = Network('MRP')
# Block 1: nodes
delay = Delay(model, 'WorkingState')
queue = Queue(model, 'RepairQueue', SchedStrategy.FCFS)
queue.setNumberOfServers(2)
# Block 2: classes
cclass = ClosedClass(model, 'Machines', 3, delay)
delay.setService(cclass, Exp(0.5))
queue.setService(cclass, Exp(4.0))
# Block 3: topology
model.link(Network.serialRouting(delay, queue))


# In[3]:


# Block 4: solution
solver = SolverCTMC(model)
ctmcAvgTable = solver.getAvgTable()
print(ctmcAvgTable)


# In[4]:


StateSpace = solver.getStateSpace()
print("\nStateSpace =")
for i in range(len(StateSpace)):
    StateSpace[i].printStateVector()


# In[5]:


InfGen = solver.getGenerator()
print("\nInfGen =")
print(InfGen)

# TBC
# model.printInfGen(InfGen,StateSpace)
#
# [StateSpace,nodeStateSpace] = solver.getStateSpace()
# nodeStateSpace{delay}
# nodeStateSpace{queue}

