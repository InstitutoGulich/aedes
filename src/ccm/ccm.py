import sys
import pylab as pl
import numpy as np
from attractor import solve
from scipy.stats import pearsonr

def exportToCsv(M,time_range):
    print('X, Y, Z')
    for i in range(0,len(time_range)):
        print('%s, %s, %s'%(M[i,0], M[i,1], M[i,2]))

#the d, and t this method returns seems to be ok, but it would be nice a second look
def nearestNeighbors(M_U,time_range,E,tau):
    d=np.zeros( (len(time_range)- (E-1)*tau,E+1) )
    idx=np.zeros( (len(time_range)- (E-1)*tau,E+1),dtype=int )
    for i in range(0,len(time_range)-(E-1)*tau):
        distances_i=np.linalg.norm(M_U-M_U[i],axis=1)
        idx[i,:]=np.argpartition(distances_i,range(E+1 +1))[1: E+1 +1]#[:E+1] -> [1: E+1 +1] to exclude distance to itself
        d[i,:]=distances_i[idx[i,:]]

    return d,idx


def createWeights(d):
    w=np.zeros( (len(time_range)- (E-1)*tau,E+1) )
    for i in range(0,len(time_range)-(E-1)*tau):
        u=np.exp(-d[i,:]/d[i,0])
        N=np.sum(u)
        w[i,:]=u/N

    return w

def estimate(M_U,V,time_range,E,tau):
    d,idx=nearestNeighbors(M_U,time_range,E,tau)
    w=createWeights(d)
    e=np.zeros( (len(time_range)- (E-1)*tau) )
    for i in range(0,len(time_range)-(E-1)*tau):
        e[i]=np.dot(w[i,:], V[idx[i,:]+(E-1)*tau])

    return e#print(e-Y[(E-1)*tau:])

def C(M_U,V,time_range,E,tau):
    e=estimate(M_U,V,time_range,E,tau)
    rho,tail=pearsonr(V[(E-1)*tau:],e)
    return rho**2

#x xmap y (y causes x)
#implementation follows: McCracken, James (2014). "Convergent cross-mapping and pairwise asymmetric inference".
if(__name__ == '__main__'):
    L_range=range(20,800,20)
    full_M,full_time_range=solve(L_range[-1],L_range[-1]*5)
    U_crossmap_V=[]
    V_crossmap_U=[]
    for L in L_range:
        idx=np.abs(full_time_range-L).argmin()
        M,time_range=full_M[:idx],full_time_range[:idx]
        tau=1
        M_X=np.array([ [M[i,0],M[i-tau,0],M[i-2*tau,0]] for i in range(2*tau,len(time_range)) ])
        M_Y=np.array([ [M[i,1],M[i-tau,1],M[i-2*tau,1]] for i in range(2*tau,len(time_range)) ])
        M_Z=np.array([ [M[i,2],M[i-tau,2],M[i-2*tau,2]] for i in range(2*tau,len(time_range)) ])
        E=3
        #exportToCsv(M,time_range)
        #nearestNeighbors(M_X,time_range,E,tau)
        #createWeights(M_X,time_range,E,tau)
        #estimate(M[:,1],M_X,time_range,E,tau)
        XYZ={'X':M[:,0],'Y':M[:,1],'Z':M[:,2]}
        M_XYZ={'X':M_X,'Y':M_Y,'Z':M_Z}
        U=XYZ[sys.argv[1].upper()]
        M_U=M_XYZ[sys.argv[1].upper()]
        V=XYZ[sys.argv[2].upper()]
        M_V=M_XYZ[sys.argv[2].upper()]
        U_crossmap_V.append(C(M_U,V,time_range,E,tau))
        V_crossmap_U.append(C(M_V,U,time_range,E,tau))
        print('.', end='', flush=True)

    pl.plot(L_range,U_crossmap_V,label='%s xmap %s'%(sys.argv[1],sys.argv[2]))
    pl.plot(L_range,V_crossmap_U,label='%s xmap %s'%(sys.argv[2],sys.argv[1]))
    pl.legend(loc=0)
    pl.show()
