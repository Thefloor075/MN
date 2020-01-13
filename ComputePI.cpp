#include <iostream>
#include <cstdlib>
#include <time.h>
#include <chrono>
#include <vector>
#include <omp.h>

#define PI 3.14159265359



int main(){
	int N = 100000;

	int i,d_N;
	int A = 10;
	//Parallele zone
	#pragma omp parallel
	{	
		#pragma omp for private(d_N,i)
		for (d_N = 10; d_N < N; d_N+=A){
			//std::chrono::steady_clock::time_point begin_t = std::chrono::steady_clock::now();		
			unsigned int s = 0;
			unsigned int tid = omp_get_thread_num();
			unsigned int tid2 = 2*tid;	
			double x,y;
			double inv = 1/RAND_MAX;
			for ( i = 0; i < d_N; i++){
				x = static_cast<double>(rand_r(&tid))*inv;
				y = static_cast<double>(rand_r(&(tid2)))*inv;
				if (x*x + y*y < 1){
					s++;
				}
			}
			if (d_N%1000 == 0){
				std::cout << 4*static_cast<double>(s)/d_N << std::endl;
			}
			//std::chrono::steady_clock::time_point end_t = std::chrono::steady_clock::now();
			//double temps = static_cast<float>(std::chrono::duration_cast<std::chrono::milliseconds>(end_t-begin_t).count())/1000;
			//std::cout << "Tps : " << temps << " N : " << d_N << std::endl;
		}
	}
}
