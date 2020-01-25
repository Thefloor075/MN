

from pylab import show,subplot,figure
import matplotlib.pyplot as plt
import numpy as np


def simule(lb,t, N0X):
	return N0X*np.exp(-lb*t)


def simule2(lb_x,lb_y,N0X,N0Y,t):
	c = N0Y - lb_x/(lb_y-lb_x)*N0X
	return c*np.exp(-lb_y*t) + N0X* lb_x/(lb_y-lb_x) * np.exp(-lb_x*t)

#Valeur initiale du nombre de noyaux de l'espèce X
depart_N0 = 100

#Nombre de simulation
Nb_simul = 50

#pas de temps
dt = 1e-3

#lambda
lb = 1/2

#Probabilité de désintégration
p = lb*dt


#Temps de simulation
temp = 6 #secondes

nb_step = int(temp/dt)

inv = 1/Nb_simul

L = np.zeros(nb_step)

vec_Time = np.linspace(0, temp, num=nb_step) #On regarde sur 0 à 6 secondes 

for _ in range(Nb_simul):
	depart = depart_N0
	for indice in range(nb_step):
		depart -= sum((np.random.rand(depart)<p).astype(int))
		L[indice] += depart



inv = 1/Nb_simul
L = inv*L

plt.title('Variation du nombre de noyaux en fonction du temps, moyenne sur {} simulation'.format(Nb_simul))
plt.grid()
plt.plot(vec_Time,simule(lb,vec_Time,depart_N0),'--',label = 'Nx(t) théorique')
plt.plot(vec_Time,L,label='Nx simulé')
plt.legend()
plt.show()


lb_x = 0.5 #1/2
lb_y = 0.6 #3/5

p_X = dt*lb_x 
p_Y = dt*lb_y

L_x = np.zeros(nb_step)
L_y = np.zeros(nb_step)

depart_N0X = 100
depart_N0Y = 50

nb_step = int(temp/dt)
for _ in range(Nb_simul):
	depart_Nx = depart_N0X
	depart_Ny = depart_N0Y
	for indice in range(nb_step):
		A = sum((np.random.rand(depart_Nx)<p_X).astype(int))
		depart_Nx -= A
		L_x[indice] += depart_Nx
		B = sum((np.random.rand(depart_Ny)<p_Y).astype(int))
		depart_Ny = depart_Ny + A - B
		L_y[indice] += depart_Ny

plt.title('Variation du nombre de noyaux en fonction du temps, moyenne sur {} simulation'.format(Nb_simul))
plt.plot(vec_Time,simule2(lb_x,lb_y,depart_N0X,depart_N0Y,vec_Time),'--', label='Ny(t) Théorique')
plt.plot(vec_Time,simule(lb_x,vec_Time,depart_N0),'--',label='Nx(t) Théorique')
plt.plot(vec_Time,L_y*inv, label='Ny(t) simulée')
plt.grid()
plt.plot(vec_Time,L_x*inv, label='Nx(t) simulée')
plt.legend()	
plt.show()






