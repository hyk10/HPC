# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 16:16:03 2016

@author: hyk10
"""
import matplotlib
import subprocess

L_d = 1.0
N_x_i = 20
T_d = 5.0
N_t_d = 50.0
alpha_d = 1.0

#subprocess.call('g++ main.cpp TriMatrix.cpp TriMatrix.h',shell=True)
#subprocess.call('./a.out',shell=True)
#subprocess.call('L_d',shell=True)
#subprocess.call('N_x_i',shell=True)
#subprocess.call('T_d',shell=True)
#subprocess.call('N_t_d',shell=True)
#subprocess.call('alpha_d',shell=True)

text_file = open("refT.txt", "r")
varTime = text_file.read().split(',')
text_file.close()

text_file = open("Nx_2.txt", "r")
varMidval = text_file.read().split(',')
text_file.close()

#format of the file is dt,dx,rmse
text_file = open("dtdx.txt", "r")
varTime = text_file.read().split(',')
text_file.close()

def plot_sol(Ue, U, x, t, st):
    fig, ax = plt.subplots(1, 1)
    for ue, u, ti in zip(Ue[st], U[st], t[st]):
        ax.plot(x, ue, 'o', label='t = {:.3f}'.format(ti))
        ax.plot(x, u)
    ax.plot(x, Ue[-1], 'o', label='t = {:.3f}'.format(t[-1]))
    ax.plot(x, U[-1])
#    handles, labels = ax.get_legend_handles_labels()
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0),
              borderaxespad=0., fontsize='medium',
              title='Numerical Solutions')  #, fancybox=True, shadow=True)   
    ax.set_title(r'Temperature distribution')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'U')