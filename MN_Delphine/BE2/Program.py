from scipy.linalg import solve
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


#Def Temps et espace
Temps = 10e-9
l = 100e-4

#Def constantes
n_x = 100
n_t = 80
dx = l/n_x
dt = Temps/n_t

w0 = 1e-3
eps_0 = 8.85e-12
eps_r = 12.6
mu_n = 0.3
mu_p = 0.015
c_n = 4.1e-14
c_p = 1.6e-14
sigma_n = 0.000784314
sigma_p = 0.000313725
e = 1.602e-19
N_d = 1e22
N_a = 5e21
N_t = 6.5e22
n_T0 = 5e21
p_T0 = 6e22
Ext = 1e6 
I_0 = 10000
e_n_th = 37.204
e_p_th = 1e-4

I_res = e_n_th * n_T0 / sigma_p / p_T0
n_0 = e_n_th*n_T0 / c_n / p_T0
p_0 = e_p_th*p_T0 / c_p / n_T0
E_0 = Ext



#Def 
E = E_0*np.ones([n_x, n_t])
n = n_0*np.ones([n_x, n_t])
p = n_0*np.ones([n_x, n_t])

#Conditions initiales
E[:,0] = Ext
n[:,0] = n_0
p[:,0] = p_0


def I(i,j):
    b = (i*dx - l/2)/w0
    return I_0 * np.exp(-b*b)


def F(Y,i,j):
    #Y[0] = E[i][j]
    #Y[1] = n[i][j]
    #Y[2] = p[i][j]
    E_ij = Y[0]
    n_ij = Y[1]
    p_ij = Y[2]
    I_ij = I(i,j)
    E_j1 = e / (eps_0 * eps_r) * ( - (mu_n * n_ij + mu_p * p_ij) * E_ij + (mu_n * n_0 + mu_p * p_0) * E_0)
    n_j1 = (e_n_th + sigma_n * I_ij)*(n_T0 - (E[i][j+1] - E[i][j-1])*0.5 / dx * eps_0 * eps_r / e) - c_n * n_ij * (N_t - (n_T0 - (E[i][j+1]-E[i][j+1]) * 0.5/dx * eps_0 * eps_r / e) ) + mu_n * n_ij * (E[i][j+1]-E[i][j-1]) * 0.5/dx + mu_n * E_ij * (n[i][j+1]-n[i][j+1]) * 0.5/dx
    p_j1 = (e_p_th + sigma_p * I_ij) * ( N_t - n_T0 - (E[i][j+1]-E[i][j-1]) * 0.5 / dx * (eps_0 * eps_r)/e) - c_p * p_ij * (n_T0 - (E[i][j+1] - E[i][j-1]) *0.5/dx * eps_0 * eps_r / e ) - mu_p * p_ij * (E[i][j+1]-E[i][j-1]) * 0.5 / dx - mu_p * E_ij * (p[i][j+1]-p[i][j-1]) * 0.5 / dx
    return np.array([E_j1, n_j1, p_j1])
	
def runge_kt4(y,i,j):
	k1 = dt*F(y,i,j)
	k2 = dt*F(y + k1*0.5,i,j)
	k3 = dt*F(y + k2*0.5,i,j)
	k4 = dt*F(y + k3,i,j)
	return y + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)

def affichage():
    Temps = np.linspace(0,(n_t-1)*dt,(n_t-1))
    Espace = np.linspace(0,(n_x)*dx,(n_x))
    Temps, Espace = np.meshgrid(Temps, Espace)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(Temps , Espace, E[:,:n_t-1], cmap=cm.coolwarm,linewidth=0, antialiased=False)
    plt.show()

def main():
    print('Running ... ')
    y0 = np.array([E_0, n_0, p_0])
    for i in range(0,n_x):
        print("    Espace = ", i)
        j = 0
        y = runge_kt4(y0,i,j)
        E[i][j] = y[0]
        n[i][j] = y[1]
        p[i][j] = y[2]
        for j in range(1,n_t-1):
            y = runge_kt4(y,i,j)
            E[i][j] = y[0]
            n[i][j] = y[1]
            p[i][j] = y[2]
    
    affichage()

    

if __name__ == '__main__':
    main()
