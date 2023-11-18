#!/usr/bin/env python
# coding: utf-8

# https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-1-a-mm1-queue

# In[8]:


from line_solver import *


# In[9]:


model = Network('M/M/1')
# Block 1: nodes
source = Source(model, 'mySource')
queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
sink = Sink(model, 'mySink')
# Block 2: classes
oclass = OpenClass(model, 'myClass')
source.setArrival(oclass, Exp(1))
queue.setService(oclass, Exp(2))
# Block 3: topology
model.link(Network.serialRouting(source, queue, sink))


# In[10]:


# Block 4: solution
AvgTable = SolverJMT(model, 'seed', 23000, 'samples', 10000).getAvgTable()


# In[11]:


# select a particular table row
ARow = tget(AvgTable, queue, oclass)
print(ARow, end='\n\n')


# In[12]:


# select a particular table row by node and class label
ARow = tget(AvgTable, 'myQueue', 'myClass')
print(ARow)

