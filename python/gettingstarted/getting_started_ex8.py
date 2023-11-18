#!/usr/bin/env python
# coding: utf-8

# https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-8-optimizing-a-performance-metric

# In[17]:


from line_solver import *
from scipy import optimize
import numpy as np


# In[18]:


def objFun(p):
    model = Network('LoadBalCQN')
    # Block 1: nodes
    delay = Delay(model, 'Think')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
    # Block 2: classes
    cclass = ClosedClass(model, 'Job1', 16, delay)
    delay.setService(cclass, Exp(1))
    queue1.setService(cclass, Exp(0.75))
    queue2.setService(cclass, Exp(0.50))
    P = model.initRoutingMatrix()
    P.set(cclass, cclass, queue1, delay, 1.0)
    P.set(cclass, cclass, queue2, delay, 1.0)
    P.set(cclass, cclass, delay, queue1, p)
    P.set(cclass, cclass, delay, queue2, 1.0 - p)
    model.link(P)
    R = SolverMVA(model, 'exact', 'verbose', False).getAvgSysRespT()
    return R


# In[39]:


p_opt = optimize.fminbound(objFun, 0, 1)
print(p_opt[0])


# In[36]:


import matplotlib.pyplot as plt
y = []
x = np.arange(0.01,1,0.01)
y = np.array(list(map(lambda x:objFun(x), x)))
plt.plot(x, y);

