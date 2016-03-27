# -*- coding: utf-8 -*-
"""
Created on Sun Mar 27 14:48:34 2016

@author: hyk10
"""
import matplotlib.pyplot as plt
import random
import math
import subprocess
import yaml 

L_d = '1.0'
theta_d = '0.5'
dt = []
dx = []
rmse = []

#compiling c++ code
for i in range(int(math.ceil(random.uniform(1.0, 1000.0)))):
    #generating random arguments but fixing L for convenience
    N_x_i = str(int(math.ceil(random.uniform(10.0, 50.0))))
    T_d = str(math.ceil(random.uniform(5.0, 50.0)))
    N_t_d = str(math.ceil(random.uniform(50.0, 10000.0)))
    alpha_d = str(math.ceil(random.uniform(1.0, 5.0)))
    #running c++ code with random numbers
    args = ['./a.out', L_d , N_x_i, T_d, N_t_d, alpha_d, theta_d]
    subprocess.call (args, shell=True)
    #importing dt,dx,rmse from c++ output txt file
    text_file = open("dtdx.txt", "r")
    varDxDt = text_file.read().split(',')
    text_file.close()
    #changing string to float format
    dt.append(yaml.load(varDxDt[0]))
    dx.append(yaml.load(varDxDt[1]))
    rmse.append(yaml.load(varDxDt[2]))

#plotting RMSE v.s dt and dx
plt.scatter(dx,rmse)
plt.title(r'RMSE v.s dx')
plt.xlabel(r'dx, unit length')
plt.ylabel(r'Root Mean Square Error')
plt.show()

plt.scatter(dt,rmse)
plt.title(r'RMSE v.s dt')
plt.xlabel(r'dt(sec)')
plt.ylabel(r'Root Mean Square Error')
plt.show()