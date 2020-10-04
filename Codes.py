
# coding: utf-8
# Wenxin Zhang
# Urban Transportation Models
# In[2]:


import numpy as np
import pandas as pd
import copy


# In[3]:


dfs = pd.read_excel('link.xlsx')
n=24
epsilon=0.00001


# In[4]:


# get demand matrix from the given string
with open('demand_rs.txt') as f:
    data0=f.read()
data=data0.split('Origin')
demand=np.zeros((n,n))
for i in range(1,n+1):
    for test in data[i].split(';'):
        try:
            j=int(test.split(':')[0].split()[-1])
            q=float(test.split(':')[1].split()[0])
            demand[i-1,j-1]=q       
        except:
            pass


# In[5]:


# 初始化
n=24 # number of nodes
graph0=np.zeros((n,n))+1000 # 图的邻接矩阵，权重是t(0)；1000表示infinity
capacity=np.zeros((n,n)) # capacity矩阵
path0=np.zeros((n,n)) # floyd算法中还原路径的辅助矩阵，tracking path


for k in range(len(dfs)): 
    i=dfs['Starting node'][k]-1 # 从0开始 0表示节点1
    j=dfs['Ending node'][k]-1
    graph0[i,j]=dfs['Free flow travel time (min)'][k]
    capacity[i,j]=dfs['Capacity (veh/hr)'][k]
    
for j in range(n):
    path0[:,j]=j


# In[6]:


# with open('test.txt') as f:
#     test=f.read()
# test1=test.split('#')


# In[7]:

# toy example
# n=13
# demand=np.zeros((n,n))
# graph0=np.zeros((n,n))+3000
# capacity=np.zeros((n,n))
# path0=np.zeros((n,n))

# for k in range(1,len(test1)):
#     i=int(re.findall(r"[-+]?\d*\.\d+|\d+",test1[k].split()[3])[0])-5
#     j=int(re.findall(r"[-+]?\d*\.\d+|\d+",test1[k].split()[4])[0])-5
#     capacity[i,j]=float(test1[k].split()[-1])
#     graph0[i,j]=float(test1[k].split()[-3][:-1])
# demand[0,10]=6000
# demand[0,12]=6750
# demand[1,10]=7500
# demand[1,12]=5250
# for j in range(n):
#     path0[:,j]=j


# In[8]:


def floyd_w(g): # Floyd-Warshall algorithm: get the shortest paths for all pairs
    path=copy.deepcopy(path0)
    graph=copy.deepcopy(g)
    for k in range(0,n):
        for i in range(0,n):
            for j in range(0,n):
                if graph[i,j] > graph[i,k] + graph[k,j]:
                    graph[i,j] = graph[i,k] + graph[k,j]
                    path[i,j] = path[i,k]
    return path


# In[9]:


# 根据floyd算法得到的path记录还原r到s的路径；将y_ij^rs加到y_ij上得到link flow
# y记录各个link flow  y_init=np.zeros((n,n))  q=q^rs demand is the matrix containing all q
def path_link(path,r,s):
    q=demand[r,s]
    y=np.zeros((n,n))
    while path[r,s]!=s:
        r_k=int(path[r,s]) # r->k->...->s
        y[r,r_k]=q
        r=r_k # keep tracking
    y[r,s]=q
    return y      


# In[10]:


def allornothing_assignment(graph):
    y=np.zeros((n,n)) # link flow
    path=floyd_w(graph)
    for i in range(n):
        for j in range(n):
            if demand[i,j]!=0:
                y=y+path_link(path,i,j)
    return y
    


# In[11]:


def update_graph(x): #note that the graph passed contain free flow travel time
    graph=copy.deepcopy(graph0)
    for i in range(n):
        for j in range(n):
            if capacity[i,j]!=0 and x[i,j]!=0:
                graph[i,j]=graph0[i,j]*(1+0.15*(x[i,j]/capacity[i,j])**4)
    return graph


# In[12]:


def bisect(func, x, y ,low=0, high=1):
    
    'Find root of continuous function where f(low) and f(high) have opposite signs'
    assert func(low,x,y)*func(high,x,y)<=0
    while(high-low>0.0001): #stopping criterion
        midpoint = (low + high) / 2.0
        if func(low,x,y)*func(midpoint,x,y)>0:
            low = midpoint
        else:
            high = midpoint

    return midpoint

def func(alpha,x,y):
    s=0
    for i in range(len(x)):
        for j in range(len(x)):
            if capacity[i,j]!=0:
                temp=(y[i,j]-x[i,j])*(graph0[i,j]*(1+0.15*((x[i,j]+alpha*(y[i,j]-x[i,j]))/capacity[i,j])**4))
                s=s+temp
    return s



# In[27]:


import time
objs=[]
time_start=time.time()
# F-W algorithm
x0=allornothing_assignment(graph0)
graph=update_graph(x0)
obj=sum(sum(x0*graph))
objs.append(obj)
y=allornothing_assignment(graph)
alpha=bisect(func, x0, y)
x1=x0+alpha*(y-x0)
count=1

while sum(sum(x1-x0)**2)/sum(sum(x0))>0.00001:
    x0=x1
    graph=update_graph(x0)
    obj=sum(sum(x0*graph))
    objs.append(obj)
    y=allornothing_assignment(graph)
    alpha=bisect(func, x0, y)
    x1=x0+alpha*(y-x0)
    count=count+1
graphf=update_graph(x1)
objs.append(sum(sum(x1*graphf)))

time_end=time.time()
print('totally cost',time_end-time_start)


# In[22]:


def update_graph_SO(x): #note that the graph passed contain free flow travel time
    graph=copy.deepcopy(graph0)
    for i in range(n):
        for j in range(n):
            if capacity[i,j]!=0 and x[i,j]!=0:
                graph[i,j]=graph0[i,j]*(1+0.75*(x[i,j]/capacity[i,j])**4)
    return graph

def func_SO(alpha,x,y):
    s=0
    for i in range(len(x)):
        for j in range(len(x)):
            if capacity[i,j]!=0:
                temp=(y[i,j]-x[i,j])*(graph0[i,j]*(1+0.75*((x[i,j]+alpha*(y[i,j]-x[i,j]))/capacity[i,j])**4))
                s=s+temp
    return s


# In[31]:


import time
objs=[]
time_start=time.time()
# F-W algorithm
x0=allornothing_assignment(graph0)
graph=update_graph_SO(x0)
g=update_graph(x0)
obj=sum(sum(x0*g))
objs.append(obj)
y=allornothing_assignment(graph)
alpha=bisect(func_SO, x0, y)
x1=x0+alpha*(y-x0)
count=1

while sum(sum(x1-x0)**2)/sum(sum(x0))>0.00001:
    x0=x1
    graph=update_graph_SO(x0)
    g=update_graph(x0)
    obj=sum(sum(x0*g))
    objs.append(obj)
    y=allornothing_assignment(graph)
    alpha=bisect(func_SO, x0, y)
    x1=x0+alpha*(y-x0)
    count=count+1


time_end=time.time()
print('totally cost',time_end-time_start)


# In[33]:


import math
import matplotlib.pyplot as plt
objso=objs[2150:]
epochs=range(1,count+2)
plt.plot(epochs,objso,color='b',label='\sum_a t_a(x_a)')
plt.xlabel('Iteration')
plt.ylabel('\sum_a t_a(x_a)')
plt.title('F-W Algorithm')
plt.legend()


# In[35]:


import math
import matplotlib.pyplot as plt
#test=[math.log(x) for x in objs]
t=objso[1500:]
epo=range(1500,2649)
plt.plot(epo,t,color='b',label='\sum_a t_a(x_a)')
plt.xlabel('Iteration')
plt.ylabel('\sum_a t_a(x_a)')
plt.title('F-W Algorithm')
plt.legend()
plt.show()


# In[49]:

# Task 2
x_so=x1
marginal=np.zeros((n,n))
for i in range(n):
    for j in range(n):
        if capacity[i,j]>0:
            marginal[i,j]=0.6*graph0[i,j]*(x_so[i,j]/capacity[i,j])**4

def update_graph_m(x): #note that the graph passed contain free flow travel time
    graph=copy.deepcopy(graph0)
    for i in range(n):
        for j in range(n):
            if capacity[i,j]!=0 and x[i,j]!=0:
                graph[i,j]=graph0[i,j]*(1+0.15*(x[i,j]/capacity[i,j])**4)+marginal[i,j]
    return graph

def func_m(alpha,x,y):
    s=0
    for i in range(len(x)):
        for j in range(len(x)):
            if capacity[i,j]!=0:
                temp=(y[i,j]-x[i,j])*(graph0[i,j]*(1+0.15*((x[i,j]+alpha*(y[i,j]-x[i,j]))/capacity[i,j])**4)+marginal[i,j])
                s=s+temp
    return s


# In[42]:

# Task 3
import time
# objs=[]
time_start=time.time()
# F-W algorithm
x0=allornothing_assignment(graph0)
graph=update_graph_m(x0)
# obj=sum(sum(x0*graph))
# objs.append(obj)
y=allornothing_assignment(graph)
alpha=bisect(func_m, x0, y)
x1=x0+alpha*(y-x0)
count=1

while sum(sum(x1-x0)**2)/sum(sum(x0))>0.00001:
    x0=x1
    graph=update_graph_m(x0)
#     obj=sum(sum(x0*graph))
#     objs.append(obj)
    y=allornothing_assignment(graph)
    alpha=bisect(func_m, x0, y)
    x1=x0+alpha*(y-x0)
    count=count+1
graphf=update_graph(x1)
obj=sum(sum(x1*graphf))

time_end=time.time()
print('totally cost',time_end-time_start,count)