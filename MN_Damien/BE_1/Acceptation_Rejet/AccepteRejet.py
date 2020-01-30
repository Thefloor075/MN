
import numpy as np
import matplotlib.pyplot as plt
from random import *

def Acceptation_rejet(n,M,C,x_min,xmax):
# M : Donné dans l'énoncé
# x_min et x_max correspond à l'intervalle sur lequel on veut simuler f
# n nombres d'échantillons
    echan = []
    for _ in range(n):
        z = random()*(x_max - x_min) + x_min
        u = np.random.uniform(0, M*g(lb,z))
        if u <= C*f(z):
            echan.append(z)

    return np.array(echan)


def f(x): 
#définition de f
	return (2 + np.sin(x)) * np.exp(- ( 2 + np.cos(3*x) + np.sin(2*x))*x)


def g(lb,x): 
#définition de g
    return lb*np.exp(-lb*x)




# intersalle sur lequel on veut simuler f 
x_min = 0
x_max = 6

#lambda
lb = 0.25

x = np.linspace(x_min,x_max,num=100000)


#supp de f(x)/g(x)
M1 = max(f(x) / g(lb,x)) 

print("M : ",M1)

#On prend M = 12
#Donné par l'énoncé
M = 12
C = 1.5


#Affichage
plt.plot(x, f(x), '--', label='f(x)')
plt.plot(x, M*g(lb,x), label = 'g(x)')
plt.xlabel('y')
plt.ylabel('x')
plt.title('Simulation de f(x) par la méthode d\'Acceptation-Rejet')
result1 = Acceptation_rejet(10000,M,C,x_min, x_max)


nb_hist = 200

plt.hist(result1, nb_hist, normed=1, facecolor='blue', alpha=0.5)
plt.legend()
plt.show()
