import random
import math
import matplotlib.pyplot as plt
import numpy as np


pi = math.pi
s = 0

A = 10
D = 10
B = 100
P = 5 



L = []

Nb_points = 1000000
rangeL = int(np.log(Nb_points))

Nb_Sim = 1000
L = [0 for _ in range(90)]
L_ind = [0 for _ in range(90)]
L_Th = [0 for _ in range(90)]
for nn in range(Nb_Sim):
	indice = 0
	s = 0
	A = 10
	D = 10
	B = 100
	P = 5 
	for i in range(D,Nb_points):
		x = random.random() 
		y = random.random()
		if (x*x + y*y <= 1):
			s+=1
		if (i == A):
			A += P
			if (A%B == 0):
				P*= 10
				B*= 10
			L_ind[indice] = i
			L_Th[indice] = i**(-0.5)
			L[indice] +=  abs(4*s/i - pi)
			indice += 1

L = (1/Nb_Sim * np.array(L)) 


plt.plot(L_ind,L)
plt.plot(L_ind,L_Th,'--')
plt.title('pi moyenné et approximation de pi en fonction du nombre de points')
plt.ylabel('Approximation du pi')
plt.xlabel('Nombre de points')
plt.legend()
plt.grid(True,which="both",ls="--")
plt.xscale('log')
plt.yscale('log')
plt.show()
