#!/usr/bin/env python
# coding: utf-8

# https://github.com/imperial-qore/line-solver/wiki/Getting-started#example-7-response-time-distribution-and-percentiles

# In[1]:


from line_solver import *
# import matplotlib.pyplot as plt


# In[2]:


model = Network('model')

# Block 1: nodes
node = np.empty(2, dtype=object)
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.PS)

# Block 2: classes
jobclass = np.empty(2, dtype=object)
jobclass[0] = ClosedClass(model, 'Class1', 5, node[0], 0)
node[0].setService(jobclass[0], Exp(1.0))
node[1].setService(jobclass[0], Exp(0.5))

# Block 3: topology
model.link(Network.serialRouting(node[0], node[1]))


# In[5]:


# Block 4: solution
RDfluid = SolverFluid(model).getCdfRespT()
print(RDfluid)


# In[6]:


RDsim = SolverJMT(model, 'seed', 23000, 'samples', 10000).getCdfRespT()
print(RDsim)

# # Plot results
# semilogx(RDsim{2,1}(:,2),1-RDsim{2,1}(:,1),'r'); hold all;
# semilogx(RDfluid{2,1}(:,2),1-RDfluid{2,1}(:,1),'k--');
# legend('jmt-transient','fluid-steady','Location','Best');
# ylabel('Pr(T > t)'); xlabel('time t');

