"""
Created on Wed Sep  2 17:05:31 2020

@author: Xinling
"""



import numpy             as np
import matplotlib.pyplot as plt
import pandas            as pd
import seaborn           as sns
from matplotlib          import rcParams
from scipy.integrate     import odeint
from scipy.optimize      import basinhopping

plt.close('all')
sns.set_context('notebook')
rcParams['figure.facecolor'] = 'white'
rcParams['axes.facecolor']   = 'white'

#model_1, only protein synthesis with inducer, no mRNA involved

def model_1step(y, t, ind, params):
    prot = y[0]
    
    inducer = ind
    k_ind   = params[0]
    synm    = params[1]
    mu      = params[2]
    
    dprot = synm*inducer/(inducer+k_ind)-mu*prot
    
    return np.array([dprot])

#Figure 1
raw_data = pd.read_csv("Example_Data.csv")
tspan    = raw_data.iloc[:,0]
y_data   = raw_data.iloc[:,1:]

fig1 = plt.figure()
ax1  = fig1.add_subplot(1,1,1)
ax1.plot(tspan, y_data,'o')
setting = ax1.set(title='Observed data', xlabel='Time (hr)', ylabel='Conc')

#initial condition
y_init  = [2.5e-8]
values  = {'k_ind'   : 12,
           'synm'    : 2e-5,
           'mu'      : 0.015
          }
params     = list(values.values())
#tspan= np.linspace(0,360)
inducer = 1

#Solve ODE and generate Fig 2
model_data = odeint(model_1step, y_init, tspan, args=tuple([inducer, params]))
y_model    = model_data[:,0]


fig2 = plt.figure()
ax2  = fig2.add_subplot(1,1,1)
ax2.plot(tspan, y_data,'o')
ax2.plot(tspan, y_model)
setting = ax1.set(title='Observed data', xlabel='Time (hr)', ylabel='Conc')

'''
#SSE
def get_SSE(tspan, y_data, inducer, params):
    model_data = odeint(model, y_init, tspan, args=tuple([inducer, params]))
    y_model    = model_data[:,1]
    
    SSE = 0
    for column in y_data:
        SSE += np.sum((y_data[column]-y_model)**2)
    
    return SSE


def objective_function_wrapper(tspan, y_data, inducer):
    def helper(params):
        return get_SSE(tspan, y_data, inducer, params)
    return helper

def take_step(params):
    delta      = np.random.rand(len(params)) -0.5
    new_params = params + delta*params
    
    return new_params

#Wrap get_SSE
objective_function = objective_function_wrapper(tspan, y_data, inducer)

result = basinhopping(objective_function, x0=params, niter=200, disp=True, 
                      take_step=take_step)
print(result)

best_value  = result.fun
best_params = list(result.x)
print(best_params)
'''
