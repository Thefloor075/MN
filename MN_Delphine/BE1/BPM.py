from scipy.linalg import solve
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


#Permet de retrouver un semblant de solution
#MAIS ce n'est pas le bon champ.

#Temps & Espace
Temps = 100e-9
X = 1e-6 * 100


#Constantes
n_X = 108
n_T = 192
dx = X/n_X
l2 = X/2
dt = 1e-9
eps_0 = 8.85e-12
eps_r = 55
nu = 6e-7
s = 3e-2
e = 1.602e-19
N_d = 3.2e21
N_a = 9e20
a_0 = 1e-5
Ext = 2500/(4e-3)
j = 0


#Matrices & Vector
E_final = Ext*np.ones([n_X, n_T])
E_0 = Ext*np.ones([n_X,1]);
E_final[:,0] = E_0[:,0]

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


def Intensity(x):
	inv = -1/(a_0*a_0)
	b = (x - l2)
	return 1e10*np.exp(inv*b*b)


def f(i,j): 
	return 1 - np.exp(- s * Intensity(dx*i) * j * dt)


def g(i,j):
	x_i = dx*i
	t_j = dt*j
	I = Intensity(x_i)
	I_1 = Intensity(x_i - dx)
	return t_j * np.exp( - s * I * t_j) * (I - I_1) * inv_dx
	

def A(i,j):
	return eps_tot_nu*E_final[i][j]

#B(i,j) est constant

def C(i,j):
	return coef_N * f(i,j)

def D(i,j):
	return coef_Ns * g(i,j)
	

def a(i,j):
	return inv_dx*(A(i,j)*inv_dx - eps_tot_invdt - C(i,j))

def b(i,j):
	return inv_dx*(- 2*A(i,j)*inv_dx + eps_tot_invdt + C(i,j)) + D(i,j)

def c(i,j):
	return A(i,j)*inv_dx*inv_dx

def d(i,j):
	#Problème sur ce coefficient. d(0,j) n'existe pas
	return eps_tot_invdt * inv_dx * (- E_final[i][j-1] + E_final[i-1][j-1])
	#_____________________________________________________________^________
	
def update_P(j):
	for k in range(n_X-1):
		P[k][k] = b(k,j)
		P[k+1][k] = a(k+1,j)
		P[k][k+1] = c(k,j)
	P[n_X-1][n_X-1] = b(n_X-1,j)

def update_B(j):
	for k in range(1,n_X-1):
		BB[k][0] = - d(k,j)
	BB[0][0] = -(a(0,j) * Ext + d(0,j))
	BB[n_X-1][0] = -(c(n_X-1,j) * Ext + d(n_X-1,j))

def update_Efinal(H,j):
	E_final[:,j] = H[:,0]


def verify_E(E):
	return np.where(E < 0, 0, E)


def affichage():
	Temps = np.linspace(0,n_T*dt,n_T)
	Espace = np.linspace(0,n_X*dx,n_X)
	Temsp, Espace = np.meshgrid(Temps, Espace)
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(Temps , Espace, E_final, cmap=cm.coolwarm,linewidth=0, antialiased=False)
	plt.show()

def main():
	error = 1e-5
	#Init
	j = 0
	update_P(j)
	update_B(j)
	for j in range(1,n_T):
		#Init Eprime
		Eprime = np.zeros([n_X,1])
		E = solve(P,BB)
		E =verify_E(E)
		E_final = update_Efinal(E,j)
		update_P(j)
		update_B(j)
		print(j)
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



	

























