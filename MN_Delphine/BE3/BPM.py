import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


#Def Temps et largeur/longueur
Temps = 100e-9
lar = 200e-5
lon = 1e-2

#Def constantes
N_Z = 300
N_lar = 300
N_T = 100
N_X = N_lar
dx = lar/N_lar
dt = Temps/N_T
dz = lon/N_Z

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

n = 3.17
e_n = 37
n_0 = e_n * n_T0 / c_n / p_T0
phi_p = (1/1.275)*1e-3

E_0 = Ext
r_eff = 1.7e-12
eps_tot = eps_0 * eps_r
A = e/eps_tot * mu_n * n_0
Q = e/eps_tot * mu_p * phi_p * p_T0 / c_p / n_T0


#longueur d'onde 
lamb_da = 1064e-9
k0 = 2*np.pi/lamb_da

#Defintion matrix and vector
F = np.zeros([N_X, N_Z])
TF_A = np.zeros([N_X, N_Z], dtype = complex)
E_final = np.zeros([N_T, N_X])
S = np.zeros([N_T, N_X], dtype=complex)
D = np.zeros([N_T, N_X], dtype=complex)



S = 1j * k0 / n


def I(i,j):
    b = (i*dx - lar/2)/w0
    return I_0 * np.exp(-b*b)


def Intensity(x):
    b = (x - lar/2)/w0
    return I_0 * np.exp(-b*b)


def operator_D():
    for i in range(1,N_T-1):
        for j in range(N_X):
            D[i,j] = 1j / (2*k0) * (E_final[i+1,j] - 2*E_final[i,j] + E_final[i-1,j]) / (dx * dx)

def operator_S():
    for i in range(N_T):
        for j in range(N_X):
            S[i,j] = 1j * k0 / n * E_final[i,j]

def E(i,j):
    In = I(i,j)
    return  E_0 / (A + Q * In) * (A + Q * In * np.exp(- (A + Q * In) * j * dt) )
    

def affichage():
    Temps = np.linspace(0,N_T*dt,N_T)
    Espace = np.linspace(0,N_X*dx,N_X)  
    Temps, Espace = np.meshgrid(Temps, Espace)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(Temps , Espace, E_final[:,:], cmap=cm.coolwarm,linewidth=0, antialiased=False)
    plt.show()


def unknow_function():
    #kx = np.linspace(-2*np.pi /lar * np.floor(N_lar/2+1),-2*np.pi /lar * np.floor(N_lar/2), N_lar)
    kx = np.linspace(-2*np.pi /lar*np.floor(N_lar/2+1),2*np.pi /lar*np.floor(N_lar/2), N_lar)
    print(kx)
    

    FFTF = np.fft.fftshift(np.exp(-1j*dz/(2*k0) * kx * kx  ))
    for z in range(N_Z-1):
        TF_A[:,z+1] = TF_A[:,z] * FFTF * np.exp(S*dz)
            

"""
for z=1:zmax-1
    FFTF=fft(F(:,z));
    FFTF=FFTF.*fftshift(exp(-i*deltaz/(2*k0)*kx.^2))';
    F(:,z+1)=ifft(FFTF).*exp(deltaz.*NL(:,z).*B(:,z));
end;
"""

def calcul_E():
    for i in range(N_T):
        for j in range(N_X):
            E_final[i,j] = E(i,j)


def affichage_2():
    Temps = np.linspace(0,N_X*dt,N_X)
    Espace = np.linspace(0,N_Z*dx,N_Z)  
    Temps, Espace = np.meshgrid(Temps, Espace)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(Temps , Espace, F[:,:], cmap=cm.coolwarm,linewidth=0, antialiased=False)
    plt.show()
    plt.plot(F[:,0])
    plt.plot(F[:,-1])
    plt.plot(F[:,0]-F[:,-1])
    plt.show()


def CI():
    Esp = np.linspace(0,N_X*dx,N_X)
    I = Intensity(Esp)
    TF_A[:,0] = np.fft.fft(np.sqrt(I))

def absF():
    for z in range(N_Z):
        F[:,z] = np.abs(np.fft.ifft(TF_A[:,z]))

def main():
    calcul_E()
    CI()
    unknow_function()
    absF()
    affichage_2()
    #unknow_function()


if __name__ == '__main__':
    main()
