from pylab import show,subplot,figure
import matplotlib.pyplot as plt
import numpy as np



def f(lb,x):
	inv = -1/lb
	return  inv * np.log(1 - x)


def simul(lb,x):
	return lb * np.exp(-lb*x)

lb = 0.5 #lambda
nb_hist	 = 100
n = 100000 #number of points
L = np.random.rand(n) #generate random numbers


F = f(lb,L)

vec_th = np.linspace(0, 15, num=n)

F_simul = simul(lb,vec_th)

plt.plot(vec_th,F_simul,'--')
plt.hist(F, nb_hist, normed=1, facecolor='red', alpha=0.5)

plt.show()


