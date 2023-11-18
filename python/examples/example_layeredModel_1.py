#!/usr/bin/env python
# coding: utf-8

# This example illustrates the execution on a layered queueing network model.
# Performance indexes now refer to processors, tasks, entries, and activities.
# Indexes refer to the submodel (layer) where the processor or task acts as a server.
# NaN indexes indicate that the metric is not supported by the node type.

# In[4]:


from line_solver import *


# In[5]:


model = LayeredNetwork('myLayeredModel')
P = np.empty(2, dtype=object)
P[0] = Processor(model, 'P1', 1, SchedStrategy.PS)
P[1] = Processor(model, 'P2', 1, SchedStrategy.PS)
T = np.empty(2, dtype=object)
T[0] = Task(model, 'T1', 10, SchedStrategy.REF).on(P[0]).setThinkTime(Exp.fitMean(100))
T[1] = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P[1]).setThinkTime(Immediate())
E = np.empty(2, dtype=object)
E[0] = Entry(model, 'E1').on(T[0])
E[1] = Entry(model, 'E2').on(T[1])
A = np.empty(4, dtype=object)
A[0] = Activity(model, 'AS1', Exp.fitMean(1.6)).on(T[0]).boundTo(E[0])
A[1] = Activity(model, 'AS2', Immediate()).on(T[0]).synchCall(E[1],1)
A[2] = Activity(model, 'AS3', Exp.fitMean(5)).on(T[1]).boundTo(E[1])
A[3] = Activity(model, 'AS4', Exp.fitMean(1)).on(T[1]).repliesTo(E[1])
# T[0].addPrecedence(ActivityPrecedence.Serial(A[0], A[1]))
# T[1].addPrecedence(ActivityPrecedence.Serial(A[2], A[3]))


# In[6]:


# options = SolverLQNS.defaultOptions
# options.keep = True # uncomment to keep the intermediate XML files generates while translating the model to LQNS
#
# solver[0] = SolverLQNS(model)
# AvgTable[0] = solver[0].getAvgTable()
# AvgTable[0]
#
# useLQNSnaming = true
# AvgTable[1] = solver[0].getAvgTable(useLQNSnaming)
# AvgTable[1]
#
#
# useLQNSnaming = true
# [AvgTable[2], CallAvgTable[2]] = solver[0].getRawAvgTables()
# AvgTable[2]
# CallAvgTable[2]
#
# AvgTable[3]=SolverLN(model).getAvgTable

