


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



# a simple Runge-Kutta integrator for multiple dependent variables and one independent variable

def rungekutta4(yprime, time, y0):
    # yprime is a list of functions, y0 is a list of initial values of y
    # time is a list of t-values at which solutions are computed
    #
    # Dependency: numpy

    N = len(time)

    y = array([thing*ones(N) for thing in y0]).T

    for ii in xrange(N-1):
        dt = time[ii+1] - time[ii]
        k1 = dt*yprime(y[ii], time[ii])
        k2 = dt*yprime(y[ii] + 0.5*k1, time[ii] + 0.5*dt)
        k3 = dt*yprime(y[ii] + 0.5*k2, time[ii] + 0.5*dt)
        k4 = dt*yprime(y[ii] + k3, time[ii+1])
        y[ii+1] = y[ii] + (k1 + 2.0*(k2 + k3) + k4)/6.0

    return y

# Miscellaneous functions
n= 1.0/3.0
kappa1 = 0.1
kappa2 = 0.1
kappa3 = 0.1
def total_energy(valpair):
    (x, y, px, py) = tuple(valpair)
    return .5*(px**2 + py**2) + (1.0/(1.0*(n+1)))*(kappa1*np.absolute(x)**(n+1)+kappa2*np.absolute(y-x)**(n+1)+kappa3*np.absolute(y)**(n+1))

def pqdot(valpair, tval):
    # input: [x, y, px, py], t
    # takes a pair of x and y values and returns \dot{p} according to the Hamiltonian
    (x, y, px, py) = tuple(valpair)
    return np.array([px, py, -kappa1*np.sign(x)*np.absolute(x)**n+kappa2*np.sign(y-x)*np.absolute(y-x)**n, kappa2*np.sign(y-x)*np.absolute(y-x)**n-kappa3*np.sign(y)*np.absolute(y)**n]).T

def findcrossings(data, data1):
    # returns indices in 1D data set where the data crossed zero. Useful for generating Poincare map at 0
    prb = list()
    for ii in xrange(len(data)-1):
        if (((data[ii] > 0) and (data[ii+1] < 0)) or ((data[ii] < 0) and (data[ii+1] > 0))) and data1[ii] > 0:
            prb.append(ii)
    return array(prb)

t = linspace(0, 1000.0, 100000)
print ("step size is " + str(t[1]-t[0]))

# Representative initial conditions for E=1
E = 1
x0=0
y0=0
init_cons = [[x0, y0, np.sqrt(2*E-(1.0*i/10.0)*(1.0*i/10.0)-2.0/(n+1)*(kappa1*np.absolute(x0)**(n+1)+kappa2*np.absolute(y0-x0)**(n+1)+kappa3*np.absolute(y0)**(n+1))), 1.0*i/10.0] for i in range(-10,11)]

outs = list()
for con in init_cons:
    outs.append( rungekutta4(pqdot, t, con) )


# plot the results
fig1 = figure(1)
for ii in xrange(4):
    subplot(2, 2, ii+1)
    plot(outs[ii][:,1],outs[ii][:,3])
    ylabel("py")
    xlabel("y")
    title("Full trajectory projected onto the plane")

fig1.suptitle('Full trajectories E = 1', fontsize=10)


# Plot Poincare sections at x=0 and px>0
fig2 = figure(2)
for ii in xrange(4):
    subplot(2, 2, ii+1)
    xcrossings = findcrossings(outs[ii][:,0], outs[ii][:,3])
    yints = [.5*(outs[ii][cross, 1] + outs[ii][cross+1, 1]) for cross in xcrossings]
    pyints = [.5*(outs[ii][cross, 3] + outs[ii][cross+1, 3]) for cross in xcrossings]
    plot(yints, pyints,'.')
    ylabel("py")
    xlabel("y")
    title("Poincare section x = 0")

fig2.suptitle('Poincare Sections E = 1', fontsize=10)

show()