#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdint.h>
#include <omp.h>
#include <math.h>

// Pour faire un programme rapide
// Ne fonctionne pas encore

#define CHANNEL_NUM 1


double Intensity(double x);

double f(int i, int j);

double g(int i,int j);

double A(int i,int j);

double C(int i,int j);

double D(int i,int j);
	
double a(int i,int j);

double b(int i,int j);

double c(int i,int j);

double d(int i,int j);

void update_P(double* P);

void update_BB(double* BB);

void update_Efinal(double* BB, int j);

void verify_E(double* E);

double norm(const double* E1,const double* E2);

void reset_E(double* E);


//Constantes
double X = pow(10.0, -6) * 100;
int j = 0;
int N_X = 80;
int N_T = 100;
double dx = X/N_X;
double l2 = X/2;
double dt = pow(10.0, -9);
double eps_0 = 8.85*pow(10.0, -12);
double eps_r = 55;
double nu = 6*pow(10.0, -7);
double s = 3*pow(10.0, -2);
double e = 1.602*pow(10.0, -19);
double N_d = 3.2*pow(10.0, 21);
double N_a = 9*pow(10.0, 20);
double a_0 = 1*pow(10.0, -5);
double Ext = 2500/(4*pow(10.0, -3));
double coef_N = e*nu* (N_d - N_a);
double coef_Ns = coef_N*s;
double inv_dx = 1/dx;
double inv_dt = 1/dt;
double eps_tot = eps_0*eps_r;
double eps_tot_nu = -eps_tot*nu;
double eps_tot_invdt = eps_tot*inv_dt;
double* E_final = (double*)malloc(N_T*N_X*sizeof(double));



int main(){

	for(int i = 0; i < N_T*N_X; i++){
		E_final[i] = Ext;	
	}

	size_t sizeVector_E = (size_t) N_X*sizeof(double); 
	double* P = (double*)malloc(N_X*N_X*sizeof(double));
	double* E = (double*)malloc(sizeVector_E);
	double* Eprime = (double*) malloc(sizeVector_E);
	double* BB = (double*) malloc(sizeVector_E);
	
	double error = pow(10.0, -5);
	double square = error*error;
	double norme;
	//Init
	int j = 0;


	update_P(P);
	update_BB(BB);
	
	std::cout << N_T << std::endl;
	for (j = 1; j < N_T; j++){

		//Init 0
		reset_E(Eprime);
		//E = solve(P,BB)
		verify_E(E);
		update_Efinal(E,j);
		update_P(P);
		update_BB(BB);
		//print(j)
		norme = norm(E,Eprime);
		while(norme > square){
			std::cout << j << " " << norme << std::endl;
			//Eprime = E;
			//E = solve(P,BB);
			verify_E(E);
			update_Efinal(E,j);
			update_P(P);
			update_BB(BB);	
			norme = norm(E,Eprime);	
		}
		update_Efinal(E,j);
	}
	std::cout << "Finish " << std::endl;
}

double Intensity(double x){
	double inv = -1/(a_0*a_0);
	double b = (x - l2);
	return pow(10.0,11.0)*exp(inv*b*b);
}

double f(int i, int j){
	return 1 - exp(- s * Intensity(static_cast<double>(dx*i)) * j * dt);
}

double g(int i,int j){
	double x_i = dx*i;
	double t_j = dt*j;
	double I = Intensity(x_i);
	double I_1 = Intensity(x_i - dx);
	return t_j * exp( - s * I * t_j) * (I - I_1) *inv_dx;
}	

double A(int i,int j){
	return eps_tot_nu*E_final[i+N_X*j];
}

double C(int i,int j){
	return coef_N * f(i,j);
}


double D(int i,int j){
	return coef_Ns * g(i,j);
}
	
double a(int i,int j){
	return inv_dx*(A(i,j)*inv_dx - eps_tot_invdt - C(i,j));
}

double b(int i,int j){
	return inv_dx*(- 2*A(i,j)*inv_dx + eps_tot_invdt + C(i,j)) + D(i,j);
}

double c(int i,int j){
	return A(i,j)*inv_dx*inv_dx;
}

double d(int i,int j){
	return eps_tot_invdt * inv_dx * (- E_final[i+N_X*(j-1)] + E_final[i-1 + N_X*(j-1)]);
}

void update_P(double* P){
	for( int k = 1; k < N_X; k++){
		P[k+N_X*k] = b(k,j);
		P[k+1 + N_X*k] = a(k,j);
		P[k+N_X*(k+1)] = c(k,j);
	}
	P[(N_X+1)*(N_X-1)] = b(N_X-1,j);
}

void update_BB(double* BB){
	for (int k = 1; k < N_X; k++){
		BB[k] = - d(k,j);
	}
	BB[0] = -(a(0,j) * Ext + d(0,j));
	BB[N_X-1] = -(c(N_X-1,j) * Ext + d(N_X-1,j));
}

void update_Efinal(double* H, int j){
	for (int i = 0; i < N_X; i++){
		E_final[i+N_X*j] = H[i];
	}	
}

void verify_E(double* E){
	for(int k = 0; k < N_X; k++){
		E[k] = std::max(E[k],0.0);
	} 
}


double norm(const double* E1,const double* E2){
	double sum = 0;
	double diff = 0;
	for (int i = 0; i < N_X; i++){
		diff = E1[i] - E2[i];
		sum += diff*diff;
	}
	return sum;
}

void reset_E(double* E){
	for(int k = 0; k < N_X; k++){
		E[k] = 0;
	}
}



