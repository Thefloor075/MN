from pylab import show,subplot,figure

import numpy as np
import matplotlib.pyplot as plt

#Transformation
def gaussian(u1,u2,std):
  pi2 = 2*np.pi #coef pour optimisation
  z1 = std*np.sqrt(-2*np.log(u1))*np.cos(pi2*u2)
  z2 = std*np.sqrt(-2*np.log(u1))*np.sin(pi2*u2)
  return z1,z2


#Simulee
def th(sigma,x):
	inv = 1/(2*sigma*sigma) #coef1 pour optimisation
	inv2 = 1/(sigma*np.sqrt((2*np.pi))) #coef2 pour optimisation
	return inv2 * np.exp( - x*x * inv)

#nombre d'échantillions
n = 10000

#création de signal aléatoires de taille n
u1 = np.random.rand(n)
u2 = np.random.rand(n)


#x_th = np.random.rand(10000)
#

x_min = -3
x_max = 3
vec_th = np.linspace(x_min, x_max, num=n)

sigma = 1

# Execution de la transformation
z1,z2 = gaussian(u1,u2,sigma)

nb_hist = 100



# Affichage


figure()
subplot(221)
plt.plot(u1)  
subplot(222)
plt.plot(u2)     
subplot(223)
plt.hist(z1, nb_hist, normed=1, facecolor='blue', alpha=0.5)
plt.plot(vec_th,th(sigma,vec_th))
subplot(224)
plt.hist(z2, nb_hist, normed=1, facecolor='green', alpha=0.5)
plt.plot(vec_th,th(sigma,vec_th))
show()

