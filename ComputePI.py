from pylab import figure


import random
import math
import matplotlib.pyplot as plt



pi = math.pi #pi
s = 0

A = 10
D = 10 #begin
B = 10
P = 5 #foot


Th = []
L = []
Err = []
piT = []
Approx_Pi = []

for i in range(D,50000000):
	#Generate 
	x = random.random() #Generate number between 0 and 1
	y = random.random()
	if (x*x + y*y <= 1): #Check if the number is in the circle
		s+=1
	#For adding point
	if (i == A):
		A += P
		if (A%B == 0):
			P*= 10
			B*= 10
		approx_pi = 4*s/i
		piT.append(pi)
		Approx_Pi.append(approx_pi)
		L.append(i) 
		Err.append(abs(pi - approx_pi)) #Relative Erro
		Th.append(i**(-0.5))
	
figure()
plt.plot(L,piT,'--')
plt.plot(L,Approx_Pi)
plt.title('pi et approximation de pi en fonction du nombre de points')
plt.xlabel('Approximation du pi')
plt.ylabel('Nombre de points')
plt.legend()
plt.grid(True,which="both",ls="--")
plt.xscale('log')
plt.yscale('log')
plt.show()

figure()
plt.title('Erreur théorique et erreur de pi en fonction du nombre de points')
plt.xlabel('Erreur')
plt.ylabel('Nombre de points')
plt.plot(L,Err, label='Erreur relative simulée')
plt.plot(L,Th,'--', label='Erreur théorique')
plt.legend()
plt.grid(True,which="both",ls="-")
plt.xscale('log')
plt.yscale('log')
plt.show()	
	


