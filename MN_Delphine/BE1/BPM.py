
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import solve
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()

from matplotlib.pyplot import imshow


#eq : e_0 e_r ( d2 E / dx dt - u E d2E / d2x ) + e u (N_d - N_A)( 1 - e^( - s I_em t ) dE/dx + e u (N_d - N_a) t * e^(-s I_em t) dIem / dx E = 0

#CI
#E(x, t) = E_ext lim x -> +infty 
#E(x, 0) = E_ext

Temps = 100e-9

n_X = 80
n_T = 100
l2 = 80e-6/2
dx = 1e-6
dt = Temps/n_T
eps_0 = 8.85e-12
eps_r = 55
nu = 6e-7
s = 3e-2
e = 1.602e-19
N_d = 3.2e21
N_a = 9e20
a_0 = 1e-5
Ext = 2500/(4e-3)


E_final = Ext*np.ones([n_X, n_T])

P = np.zeros([n_X,n_X])

E = np.zeros([n_X,1])
Eprime = np.zeros([n_X,1])

BB = np.zeros([n_X,1])

coef_N = e*nu* (N_d - N_a)
coef_Ns = coef_N*s

inv_dx = 1/dx
inv_dt = 1/dt

eps_tot = eps_0*eps_r
eps_tot_nu = -eps_tot*nu

def Intensity(x):
	inv = -1/(a_0*a_0)
	b = (x - l2)
	return 1e10*np.exp(inv*b*b)

def f(i,j):
	I = Intensity(dx*i) 
	return 1 - np.exp(- s * I * j * dt)

def g(i,j):
	x_i = dx*i
	t_j = dt*j
	#f(x) = e^ax df/dx = a e^ax
	I = Intensity(x_i)
	I_1 = Intensity(x_i - dx)
	return t_j * np.exp( - s * I * t_j) * (I - I_1) / dx
	

def A(i,j):
	return eps_tot_nu*E_final[i][j]

def C(i,j):
	return coef_N * f(i,j)

def D(i,j):
	return coef_Ns * g(i,j)
	
	
def a(i,j):
	return inv_dx*(A(i,j)*inv_dx - eps_tot*inv_dt - C(i,j))

def b(i,j):
	return inv_dx*(- 2*A(i,j)*inv_dx + eps_tot*inv_dt + C(i,j)) + D(i,j)

def c(i,j):
	return A(i,j)*inv_dx*inv_dx

def d(i,j):
	return eps_tot * inv_dx * inv_dt * (- E_final[i][j-1] + E_final[i-1][j-1])

def update_P(P):
	for k in range(n_X-1):
		P[k][k] = b(k,j)
		P[k+1][k] = a(k,j)
		P[k][k+1] = c(k,j)
	P[n_X-1][n_X-1] = b(n_X-1,j)

def update_B(BB):
	#k for X
	for k in range(1,n_X-1):
		BB[k][0] = - d(k,j)
	BB[0][0] = -(a(0,j) * Ext + d(0,j))
	BB[n_X-1][0] = -(c(n_X-1,j) * Ext + d(n_X-1,j)) 

def update_Efinal(H,j):
	E_final[:,j] = H[:,0]


def verify_E(E):
	for i in range(n_X):
		if E[i] < 0:
			E[i] = 0 

error = 1e-8

#Init
j = 0


update_P(P)
update_B(BB)


for j in range(1,n_T):
	#Init 0
	Eprime = np.zeros([n_X,1]) 
	E = solve(P,BB)
	verify_E(E)
	update_Efinal(E,j)
	update_P(P)
	update_B(BB)	
	while np.linalg.norm(Eprime - E) > error:
		Eprime = E
		E = solve(P,BB)
		verify_E(E)
		update_Efinal(E,j)
		update_P(P)
		update_B(BB)	
	update_Efinal(E,j)

imshow(E_final)
plt.show()

	












