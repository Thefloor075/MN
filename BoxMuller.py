from pylab import show,subplot,figure

import numpy as np
import matplotlib.pyplot as plt

# transformation function
def gaussian(u1,u2,std):
  pi2 = 2*np.pi
  z1 = std*np.sqrt(-2*np.log(u1))*np.cos(pi2*u2)
  z2 = std*np.sqrt(-2*np.log(u1))*np.sin(pi2*u2)
  return z1,z2



def th(sigma,x):
	inv = 1/(2*sigma*sigma)
	inv2 = 1/(sigma*np.sqrt((2*np.pi)))
	return inv2 * np.exp( - x*x * inv)


# uniformly distributed values between 0 and 1
u1 = np.random.rand(10000)
u2 = np.random.rand(10000)


#x_th = np.random.rand(10000)
vec_th = np.linspace(-3, 3, num=10000)

sigma = 1

# run the transformation
z1,z2 = gaussian(u1,u2,sigma)

nb_hist = 100

#Theory
#n, bins, patches = plt.hist(z1, nb_hist, normed=1, facecolor='blue', alpha=0.5)
#show()

# plotting the values before and after the transformation




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

