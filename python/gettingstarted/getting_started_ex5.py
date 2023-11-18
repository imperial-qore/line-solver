#!/usr/bin/env python
# coding: utf-8

# https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-5-modelling-a-re-entrant-line

# In[14]:


from line_solver import *


# In[15]:


model = Network('RL')
queue = Queue(model, 'Queue', SchedStrategy.PS)
K = 3
N = (1, 0, 0)
jobclass = []
for k in range(K):
    jobclass.append(ClosedClass(model, 'Class' + str(k), N[k], queue))
    queue.setService(jobclass[k], Erlang.fitMeanAndOrder(k, 2))
P = model.initRoutingMatrix()
P.set(jobclass[0], jobclass[1], queue, queue, 1.0)
P.set(jobclass[1], jobclass[2], queue, queue, 1.0)
P.set(jobclass[2], jobclass[0], queue, queue, 1.0)
model.link(P)


# In[16]:


ctmcAvgTable = SolverMVA(model).getAvgTable()
print(ctmcAvgTable, end="\n\n")


# In[17]:


ctmcAvgSysTable = SolverMVA(model).getAvgSysTable()
print(ctmcAvgSysTable, end="\n\n")


# In[18]:


jobclass[0].completes = False
jobclass[1].completes = False
ctmcAvgSysTable2 = SolverMVA(model).getAvgSysTable()
print(ctmcAvgSysTable2)

