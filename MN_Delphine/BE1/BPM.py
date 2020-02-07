from scipy.linalg import solve
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

#Temps & Espace
Temps = 50e-9
largeur = 80e-6

#Constantes
n_X = 80
n_T = 100
dx = largeur/n_X
l2 = largeur/2
dt = Temps/n_T
eps_0 = 8.85e-12
eps_r = 56
nu = 1e-5
s = 2e-5
e = 1.602e-19
N_d = 1e25
N_a = 1e22
a_0 = 1e-5
Ext = 2500/(4e-3)
I0 = 1e10
j = 0


#Matrices & Vector
E_final = Ext*np.ones([n_X, n_T])

P = np.zeros([n_X,n_X])
E = np.zeros([n_X,1])

Eprime = np.zeros([n_X,1])
BB = np.zeros([n_X,1])

#Coefficients pour optimisation
coef_N = e*nu* (N_d - N_a)
coef_Ns = coef_N*s

inv_dx = 1/dx
inv_dt = 1/dt

eps_tot = eps_0*eps_r
eps_tot_nu = -eps_tot*nu

eps_tot_invdt = eps_tot*inv_dt


def Intensity(i,j):
	b = (i*dx - l2)/a_0
	return I0*np.exp(-b*b)

def d_Intensity(i,j):
	return -2*(i*dx-l2)/(a_0**2) * Intensity(i,j)

def f(i,j): 
	return 1 - np.exp(- s * Intensity(i,j) * j * dt)

def g(i,j):
	t_j = j*dt
	return t_j * d_Intensity(i,j) * np.exp( - s * Intensity(i,j) * t_j ) 
	
def A(i,j):
	return eps_tot_nu*E_final[i][j]

def B(i,j):
	return eps_tot

def C(i,j):
	return  e * nu * (N_d - N_a) * f(i,j)

def D(i,j):
	return  e * nu * (N_d - N_a) * s * g(i,j)
	

def a(i,j):
	return inv_dx * ( A(i,j) * inv_dx - B(i,j) * inv_dt - C(i,j))

def b(i,j):
	return inv_dx*(- 2*A(i,j) * inv_dx + B(i,j) * inv_dt + C(i,j)) + D(i,j)

def c(i,j):
	return A(i,j) * inv_dx * inv_dx

def d(i,j):
	if i == 0:
		return B(i,j) * inv_dx * inv_dt * (- E_final[i][j-1] + Ext)
	return  B(i,j) * inv_dx * inv_dt * (- E_final[i][j-1] + E_final[i-1][j-1])

	
def update_P(j):
	for k in range(0,n_X-1):
		P[k][k] = b(k,j)
		P[k+1][k] = a(k+1,j)
		P[k][k+1] = c(k,j)
	P[n_X-1][n_X-1] = b(n_X-1,j)

def update_B(j):
	for k in range(1,n_X-1):
		BB[k][0] = - d(k,j)
	BB[0][0] = -(a(0,j) * Ext + d(0,j))
	BB[n_X-1][0] = -(c(n_X-1,j) * Ext + d(n_X-1,j))

def update_Efinal(E,j):
	E_final[:,j] = E[:,0]

def verify_E(E):
	return np.where(E < 0, 0, E)

def affichage():
	Temps = np.linspace(0,n_T*dt,n_T)
	Espace = np.linspace(0,n_X*dx,n_X)
	Temps, Espace = np.meshgrid(Temps, Espace)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(Temps , Espace, E_final, cmap=cm.coolwarm,linewidth=0, antialiased=False)
	plt.show()

def main():
	#Critère de convergence
	error = 1e-7

	#On commence à T = 1, T = 0 est la condition initiale
	for j in range(1,n_T):
		#Init Eprime

		Eprime = np.zeros([n_X,1])

		#On actualise P et le vecteur B
		update_P(j)
		update_B(j)
				
		E = solve(P,BB)
		#Dans le cas ou des valeurs sont inférieurs à 0, ce qui n'est pas physique.
		E = verify_E(E)

		#On écrit le champ qu'on à calculer dans le solution final.
		E_final = update_Efinal(E,j)
		
		#On recalcul la matrice P et le vector B
		update_P(j)
		update_B(j)
		print('T = ', round(j*dt,10))
		while np.linalg.norm(Eprime - E) > error:
			Eprime = E
			E = solve(P,BB)
			E = verify_E(E)
			E_final = update_Efinal(E,j)
			update_P(j)
			update_B(j)	
		update_Efinal(E,j)
	affichage()
	
if __name__ == "__main__":
	main()



	

























