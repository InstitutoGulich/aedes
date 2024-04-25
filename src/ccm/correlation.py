import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import stats

time_range=np.array(range(0,200))
X=np.array([0.]* len(time_range))
Y=np.array([0.]* len(time_range))

X[0]=0.2
Y[0]=0.4
#https://www.nature.com/articles/srep14750
#http://www.science.ntu.edu.tw/upload/news/causalityCCM.pdf
for t in time_range[:-1]:
    X[t+1]=X[t] * (3.8-3.8*X[t]-0.02*Y[t])
    Y[t+1]=Y[t] * (3.5-3.5*Y[t]-0.1*X[t])

window=range(9,17)#(9,30)#100,130
plt.plot(time_range[window],X[window],'-o',label='X(t)')
plt.plot(time_range[window],Y[window],'-o',label='Y(t)')
plt.legend(loc=0)
print('(Pearson s correlation coefficient,2-tailed p-value): ' + str(stats.pearsonr(X[window],Y[window])) )
plt.show()
