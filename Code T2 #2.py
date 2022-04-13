


# Course: Métodos Computacionais da Física B (FIS01206)


## Note: I decided to write the paper in English as training. I affirm that there may be grammatical errors and I kindly request that they be disregarded.


# Instructions:
'''
Reproduzam os gráficos da página 167 do livro Classical Dynamics of Particles and Systems, (5a ed, aut. Marion, Thornton). 

Use Runge-Kutta 4ª ordem clássico para resolver a equação

    d^2x/dt^2 = -c dx/dt - sin (x) + F cos (wt) 

onde \quad x(0) = 1.0, v(0) = 0.0. 

Escreva a equação de 2a ordem como o sistema de equações de 1ª ordem.

    dx/dt = v 
	dv/dt = -c v - sin (x) + F cos (wt)

com  h=0.01  (x: posição angular, v: velocidade angular). Os outros parâmetros são: c=0.05, w=0.7.
'''

#___________________________________________________________

# Imports

from telnetlib import AYT
import matplotlib.pyplot as plt
import math, operator, random
from numpy import pi, sin, cos, linspace
from copy import deepcopy

#___________________________________________________________

# Runge-Kutta 4 method:

dataX0 = []
dataY0 = []
vList = []
xList = [] 

count = 1
t = 0

x = 1.0
v = 0
f = 0
z = 0

xList.append(x)
vList.append(v)

h = 0.01

c = 0.05
w = 0.7


def dxdt(x):
	dxdt = v
	return dxdt



def dvdt(x):
	dvdt = -c*v - sin(x) + f*cos(w*t)
	return dydt



def rk4(x, y):

	global h

	k1x = h*dxdt(x)
	k1y = h*dvdt(x)
	
	k2x = h*dxdt(x + k1x/2.0, y + k1y/2.0)
	k2y = h*dvdt(x + k1x/2.0, y + k1y/2.0)
	
	k3x = h*dxdt(x + k2x/2.0, y + k2y/2.0)
	k3y = h*dvdt(x + k2x/2.0, y + k2y/2.0)
	
	k4x = h*dxdt(x + k3x, y + k3y)
	k4y = h*dvdt(x + k3x, y + k3y)
	
	x = x + k1x/6.0 + k2x/3.0 + k3x/3.0 + k4x/6.0
	v = v + k1y/6.0 + k2y/3.0 + k3y/3.0 + k4y/6.0

	return [x,v]

X0 = []
V0 = []