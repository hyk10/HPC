# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 16:16:03 2016

@author: hyk10
"""
import matplotlib.pyplot as plt
import subprocess
import yaml 
import time

L_d = '1.0'
N_x_i = '20'
T_d = '5.0'
N_t_d = '5000.0'
alpha_d = '1.0'

#complie and running c++ code
args = ['g++ main.cpp TriMatrix.cpp TriMatrix.h', './a.out', L_d , N_x_i, T_d, N_t_d, alpha_d]
subprocess.call (args, shell=True)
print('hello')
time.sleep(50)
print('awake')
#importing time from c++ output txt file
text_file = open("refT.txt", "r")
varTime = text_file.read().split(',')
text_file.close()
varTime_n = []
#changing string to float format
for x in varTime:
    varTime_n.append(yaml.load(x)) 

#importing midpoint value from c++ output txt file
text_file = open("Nx_2.txt", "r")
varMidval = text_file.read().split(',')
text_file.close()
varMidval_n = []
#changing string to float format
for x in varMidval:
    varMidval_n.append(yaml.load(x)) 

#plotting Mid Value($U_{0.5L}$) v.s. Time
plt.plot(varTime_n,varMidval_n)
plt.title(r'Mid Value($U_{0.5L}$) v.s. Time')
plt.ylabel(r'Mid Value($U_{0.5L}$)')
plt.xlabel(r'Time(sec)')
plt.show()