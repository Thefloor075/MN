from scipy.linalg import solve
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


#Def Temps et largeur/Longueur
Temps = 100e-9
lar = 80e-6
lon = 1

#Def constantes
N_lar = 100
N_T = 80
N_X = 100
dx = lar/N_X
dt = Temps/N_T

w0 = 2e-5
eps_0 = 8.85e-12
eps_r = 12.6
mu_n = 0.4
mu_p = 0.015
e_n_ = 37
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

n_0 = 3.17
phi_p = (1/1.275)*1e-3

E_0 = Ext
r_eff = 1.7e-12
eps_tot = eps_0 * eps_r
A = e/eps_tot * mu_n * n_0
Q = e/eps_tot * mu_p * phi_p * p_T0 / c_p / n_T0


#Def matrices

E_final = np.zeros([N_X, N_T])

#k_x = np.linspace(-2*np.pi /(lar/2), -2*np.pi /(lar/2), N_lar)

def I(i,j):
    b = (i*dx - lar/2)/w0
    return I_0 * np.exp(-b*b)

def E(i,j):
    In = I(i,j)
    B = E_0 / (A + Q * In) * (A + Q * In * np.exp(- (A + Q * In) * j * dt) )
    return B

def affichage():
    Temps = np.linspace(0,N_T*dt,N_T)
    Espace = np.linspace(0,N_X*dx,N_X)  
    Temps, Espace = np.meshgrid(Temps, Espace)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(Temps , Espace, E_final[:,:], cmap=cm.coolwarm,linewidth=0, antialiased=False)
    plt.show()

def main():
    for i in range(0,N_X):
        for j in range(0, N_T):
            E_final[i][j] = E(i,j)
    print(E_final)
    affichage()


if __name__ == '__main__':
    main()