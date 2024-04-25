#coding: utf-8
import sys
import numpy as np
from math import sin
import scipy.integrate as spi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
def lorenz(V,t):#[-8,8,27] parameters taken from http://node99.org/tutorials/ar/
    delta=10.
    beta=8/3.
    rho=28.
    x,y,z=V
    return [delta*(y-x),x*(rho-z)-y,x*y-beta*z]

def halvorsen(V,t):#[0,1,0]
    a=1.89
    x,y,z=V
    return [-a*x -4*y -4*z - y**2, -a*y - 4*z - 4*x -z**2,-a*z - 4*x - 4*y - x**2]

def thomas(V,t):#[0,0,0.5]
    b=0.208186
    x,y,z=V
    return [sin(y)-b*x,sin(z)-b*y,sin(x)-b*z]

def logistic(V,t):#[0.2,0.4,0]
    x,y,z=V
    return [x * (3.8-3.8*x-0.02*y), y * (3.5-3.5*y-0.1*x),0]

def diff_eqs(V,t):
    '''The main set of equations'''
    return lorenz(V,t)

def solve(L=100,num=100**2):
    initial_condition=[-8,8,27]
    time_range=np.linspace(0,L,num)
    return spi.odeint(diff_eqs,initial_condition,time_range),time_range


def getCurves(fig,M,M_X,M_Y,M_Z):
    ax=Axes3D(fig)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    curve_M, = ax.plot(M[:,0], M[:,1],M[:,2],'-o',markevery=[-1],label='= (X(t), Y(t), Z(t))')
    curve_M_X, = ax.plot(M_X[:,0], M_X[:,1],M_X[:,2],'-o',markevery=[-1],label=r'= (X(t), X(t-$\tau$),X(t-2$\tau$))')
    curve_M_Y, = ax.plot(M_Y[:,0], M_Y[:,1],M_Y[:,2],'-o',markevery=[-1],label=r'= (Y(t), Y(t-$\tau$),Y(t-2$\tau$))')
    curve_M_Z, = ax.plot(M_Z[:,0], M_Z[:,1],M_Z[:,2],'-o',markevery=[-1],label=r'= (Z(t), Z(t-$\tau$),Z(t-2$\tau$))')
    line,=ax.plot([],[])
    return curve_M,curve_M_X,curve_M_Y,curve_M_Z,line

def animate(i,M,M_X,M_Y,M_Z,curve_M,curve_M_X,curve_M_Y,curve_M_Z,line):
    curve_M.set_data(M[:i,0],M[:i,1])  # update the data
    curve_M.set_3d_properties(M[:i,2])
    curve_M_X.set_data(M_X[:i,0],M_X[:i,1])  # update the data
    curve_M_X.set_3d_properties(M_X[:i,2])
    curve_M_Y.set_data(M_Y[:i,0],M_Y[:i,1])  # update the data
    curve_M_Y.set_3d_properties(M_Y[:i,2])
    curve_M_Z.set_data(M_Z[:i,0],M_Z[:i,1])  # update the data
    curve_M_Z.set_3d_properties(M_Z[:i,2])
    line.set_data([ M[i-1,0],M_X[i-1,0] ],[ M[i-1,1],M_X[i-1,1] ])
    line.set_3d_properties([ M[i-1,2],M_X[i-1,2] ])


if(__name__ == '__main__'):
    fig = plt.figure()
    M,time_range=solve()
    #Animation code based on https://matplotlib.org/examples/animation/simple_anim.html
    tau=2
    M_X=np.array([ [M[i,0],M[i-tau,0],M[i-2*tau,0]] for i in range(2*tau,len(time_range)) ])#+ np.array([20,0,0])
    M_Y=np.array([ [M[i,1],M[i-tau,1],M[i-2*tau,1]] for i in range(2*tau,len(time_range)) ])#+ np.array([-20,0,0])
    M_Z=np.array([ [M[i,2],M[i-tau,2],M[i-2*tau,2]] for i in range(2*tau,len(time_range)) ])#+ np.array([-20,0,0])
    curve_M,curve_M_X,curve_M_Y,curve_M_Z,line=getCurves(fig,M,M_X,M_Y,M_Z)
    ani = animation.FuncAnimation(fig, animate, range(0,len(time_range)), interval=100,fargs=(M,M_X,M_Y,M_Z,curve_M,curve_M_X,curve_M_Y,curve_M_Z,line))

    #deal with parameters
    curves={'M_X':curve_M_X,'M_Y':curve_M_Y,'M_Z':curve_M_Z}

    for curve_name in curves.keys():
        if(curve_name not in ['M_'+variable_name.upper() for variable_name in sys.argv[1:] ] ):
            curves[curve_name].set_visible(False)
            curves[curve_name].set_label(None)

    if('M_X' not in ['M_'+variable_name.upper() for variable_name in sys.argv[1:] ]): line.set_visible(False)#TODO:check how to make this configurable by input
    plt.legend(loc=0)
    plt.show()
