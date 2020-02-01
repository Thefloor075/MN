#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdint.h>
#include <omp.h>
#include <math.h>
#include <Eigen/Core> 
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <bits/stdc++.h>

// Pour faire un programme rapide
// Ne fonctionne pas encore

using namespace Eigen;

int N_X;
int N_T;
int nb_thread;

int array_size(int &N_X, int &N_T, int &nb_thread)
{
	std::ifstream input("input.txt");
	input >> N_X >> N_T >> nb_thread;
	return 1;
	
}
int* array = new int[array_size(N_X, N_T, nb_thread)];

// Vector & Matrix

MatrixXd E_final(N_X,N_T);
MatrixXd P(N_X,N_X);
VectorXd E(N_X);
VectorXd Eprime(N_X);
VectorXd BB(N_X);

//Constantes

double X = pow(10.0, -6) * 100;
int j = 1;

double t = pow(10.0, -9)* 100;

double dx = X/N_X;
double l2 = X/2;
double dt = t/N_T;
double eps_0 = 8.85*pow(10.0, -12);
double eps_r = 55;
double nu = 6*pow(10.0, -7);
double s = 3*pow(10.0, -2);
double e = 1.602*pow(10.0, -19);
double N_d = 3.2*pow(10.0, 21);
double N_a = 9*pow(10.0, 20);
double a_0 = 1*pow(10.0, -5);
double Ext = 625000;

double coef_N = e*nu* (N_d - N_a);
double coef_Ns = coef_N*s;
double inv_dx = 1/dx;
double inv_dt = 1/dt;
double eps_tot = eps_0*eps_r;
double eps_tot_nu = -eps_tot*nu;
double eps_tot_invdt = eps_tot*inv_dt;

inline double Intensity(double x);

inline double f(int i, int j);

inline double g(int i,int j);

inline double A(int i,int j);

inline double C(int i,int j);

inline double D(int i,int j);
	
inline double a(int i,int j);

inline double b(int i,int j);

inline double c(int i,int j);

inline double d(int i,int j);

inline void update_P_BB(MatrixXd& P, VectorXd& BB);

inline void update_Efinal(VectorXd& BB, int j);

inline void reset_E(VectorXd& E);

int main(){
	delete array;

	std::ofstream ppmFile("output_image.ppm", std::ios::out | std::ios::binary);
	ppmFile << "P6" << std::endl;
	ppmFile << N_T << " " << N_X << std::endl;
	ppmFile << "255" << std::endl;

	omp_set_num_threads(nb_thread);
	Eigen::setNbThreads(nb_thread);

	const int long WIDTH = N_T;
	const int long HEIGHT = N_X;

	#pragma omp parallel for
	for (int p = 0; p < N_X; p++){
		E(p) = 0;
		Eprime(p) = 0;
		BB(p) = 0;
		for(int c = 0; c < N_X; c++){
			P(p,c) = 0;	
		}
	}
	
	#pragma omp parallel for
	for (int p = 0; p < N_X; p++){
		for(int c = 0; c < N_T; c++){
			E_final(p,c) = Ext;			
		}
	}
	double error = pow(10.0, -3);
	//Init

	update_P_BB(P,BB);

	for (; j < N_T; j++){
		reset_E(Eprime);
		E = P.inverse()*BB;
		update_Efinal(E,j);
		update_P_BB(P,BB);
		std::cout << "j : " << j << std::endl;
		while((E-Eprime).norm() > error){
			Eprime = E;
			E = P.inverse()*BB;
			update_Efinal(E,j);
			update_P_BB(P,BB);
			
		}
	}
	
	
	int min = Ext;
	for (int x = 0; x < HEIGHT; x++){
		for (int y = 0; y < WIDTH; y++){
			if(E_final(x,y) < min){
				min = E_final(x,y);
			}
		}	
	}
	uint8_t rgbTriplet[3];
	for (int x = 0; x < HEIGHT; x++){
		for (int y = 0; y < WIDTH; y++){
			    int n = E_final(x,y)/static_cast<int>(E_final(x,y)-min);
			    rgbTriplet[0] = sqrt(n);
			    rgbTriplet[1] = n;
			    rgbTriplet[2] = n;
			    ppmFile.write((char *)rgbTriplet, 3);
		}

	}
	//stbi_write_jpg("stbjpg3.jpg", WIDTH, HEIGHT, CHANNEL_NUM, rgb_image, WIDTH);
	//free(rgb_image);
	return 0;


}

inline double Intensity(double x){
	static double I = pow(10.0,11.0);
	static double inv = -1/(a_0*a_0);
	double b = (x - l2);
	return I*exp(inv*b*b);
}

inline double f(int i, int j){
	return 1 - exp(- s * Intensity(dx*i) * j * dt);
}

inline double g(int i,int j){
	double x_i = dx*i;
	double t_j = dt*j;
	double I = Intensity(dx*i);
	return t_j * exp( - s * I * t_j) * (I - Intensity(x_i - dx)) *inv_dx;
}	

inline double A(int i,int j){
	return eps_tot_nu*E_final(i,j);
}

inline double C(int i,int j){
	return coef_N * f(i,j);
}


inline double D(int i,int j){
	return coef_Ns * g(i,j);
}
	
inline double a(int i,int j){
	return inv_dx*(A(i,j)*inv_dx - eps_tot_invdt - C(i,j));
}

inline double b(int i,int j){
	return inv_dx*(- 2*A(i,j)*inv_dx + eps_tot_invdt + C(i,j)) + D(i,j);
}

inline double c(int i,int j){
	return A(i,j)*inv_dx*inv_dx;
}

inline double d(int i,int j){
	if( i == 0){
		return eps_tot_invdt * inv_dx * (- E_final(i,j-1) + Ext);
	}
	return eps_tot_invdt * inv_dx * (- E_final(i,j-1) + E_final(i-1,j-1));
}

inline void update_P_BB(MatrixXd& P, VectorXd& BB){
	#pragma omp parallel for
	for( int k = 1; k < N_X-1; k++){
		BB(k) = - d(k,j);
		P(k,k) = b(k,j);
		P(k+1,k) = a(k+1,j);
		P(k,k+1) = c(k,j);
	}
	P(0,0) = b(0,j);
	P(1,0) = a(1,j);
	P(0,0+1) = c(0,j);
	P(N_X-1,N_X-1) = b(N_X-1,j);
	BB(0) = -(a(0,j) * Ext + d(0,j));
	BB(N_X-1) = -(c(N_X-1,j) * Ext + d(N_X-1,j));
}


inline void update_Efinal(VectorXd& H,int j){
	for (int i = 0; i < N_X; i++){
		double A = H(i);
		H(i) = std::max(H(i), 0.0);
		E_final(i,j) = A;
	}	
}

inline void reset_E(VectorXd& E){
	for(int i = 0; i < N_X; i++){
		E(i) = 0;
	}
}
