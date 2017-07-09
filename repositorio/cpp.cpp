//1	Código C++ para realizar cálculos
//1.1	Encabezado
#include <iostream>
#include <ctime>
#include <time.h>  
#include <random>
#include <fstream>
#include <iomanip>
#include <new>
#include <math.h>
#include <thread>
#include <mutex>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <string>
#include <algorithm>
#include "klu.h"

using namespace std;

/* Parámetros por defecto*/
const int N = 1000; //número de nanotubos; valor de referencia:10^3 - 10^4
const double L = 1000; //largo nanotubo [nm]
const int n = 9; //numero de segmento del nanotubo
const int n_aglomerados = 0; // numero de aglomeraos aleatorios; 0 para aglomerados fijos
const int n_aglomerados_x = 2; //numero de aglomerados fijos en x;
const int n_aglomerados_y = 1; //numero de aglomerados fijos en y;
const int n_aglomerados_z = 1; //numero de aglomerados fijos en z;
const double p_aglomerado = 0; //probabilidad de que un nanotubo este en un aglomerado
const double pi = std::acos(-1.0);
const double theta_max = (pi/2); //rango es (-theta_max,theta_max); pi/2 aleatorio; 0 alineado
const double d0 = 25; //diametro nanotubo [nm]
const double dvdw = 1.4; //distancia de equilibrio de van der waals [nm]
const double dt = 4; //distancia cut-off del efecto tunel
const double lambda =  0.2; //altura de la barrera de potencial del efecto tunel
const double R_pol = 10001; //resistencia por defecto en caso de no ocurrir percolacion

const double Lx = 5*L; //largo x del volumen representativo; valor referencia: 5
const double Ly = 5*L; //largo y del volumen representativo; valor referencia: 5
const double Lz = 5*L; //largo z del volumen representativo; valor referencia: 5

const double strain = 0;
const double poisson  = 0.45;
const int N_MonteCarlo = 1000; //numero de veces que se repite el calculo: velor de referencia: 10^2 - 10^3
const int first_seed = 1; //valor de la primera semilla
const int guardar = 2;//1 si se guardan todas las variables,2 si se guardan los resumenes, 0 si no se guarda nada
const int imprimir = 2;//1 si se guardan imprime el resumen de los resultados, 2 si se imprime los resultados de cada semilla, 0 si no se imprimen
const int resistencias = 1 ;//1 si se calcula la resistencia, 0 si no
const int porcentaje_contactos = 5; //porcentaje máximo de contactos entre nanotubos
std::mutex myMutex;

/*header de funciones implementadas*/
int redondear(double x);
double rand_xor128(int seed,long llamadas_por_seed); //generador de numeros aleatorios en (0,1). En Matlab se implemetó el mismo generador
double dotp(double a[],double b[],int size); //funcion que calcula el producto escalar entre dos array
double mean(int a[],int size); //funcion que calcula el promedio de los elementos de un array tipo int
double mean(double a[],int size); //funcion que calcula el promedio de los elementos de un array tipo double
double sd(int a[],int size); //funcion que calcula el promedio de los elementos de un array tipo int
double sd(double a[],int size); //funcion que calcula el promedio de los elementos de un array tipo double
void imprimir1d(double myArray[], int size); //funcion que imprime en pantalla un array de una dimension
void imprimir1d(int myArray[], int size); //funcion que imprime en pantalla un array de una dimension
void imprimir2d(double matrix[], int rows, int cols); //funcion que imprime en pantalla un array de dos dimensiones
void imprimir2d(int matrix[], int rows, int cols); //funcion que imprime en pantalla un array de dos dimensiones
void imprimir2dmin(double matrix[], int rows, int cols); //funcion que imprime en pantalla un array de dos dimensiones
void guardar1d(std::ofstream &outputFile,int myArray[], int size); //funcion que guarda en archivo un array tipo int de una dimension
void guardar1d(std::ofstream &outputFile,double myArray[], int size); //funcion que guarda en archivo un array tipo double de una dimension
void guardar2d(std::ofstream &outputFile,double matrix[], int rows, int cols); //funcion que guarda en archivo un array de dos dimensiones
void guardar_constantes(std::ofstream &outputFile,int seed,int N,double L,int n, int n_aglomerados,int n_aglomerados_x,int n_aglomerados_y,int n_aglomerados_z,double p_aglomerado,double theta_max,double d0,double dt, double Lx,double Ly, double Lz,double strain, double poisson);

double G_tunel(double d); //funcion que calcula la conductancia por efecto tunel a partir de la distancia entre dos cnt
void bubble_sort(int Ai[], double Ax[], int ni, int nf);
void bubble_sortA(double A_suma[], int numero_nodos);
void spice(std::ofstream &outputFile,double dmin_array[], int length_bb,int nodo_extremo[],int N,double lambda); //funcion que crea archivo ingresado a SPICE

void configuracion_inicial(double x[], double y[],double z[],double theta_suma[],double phi_suma[],double xcm[],double ycm[],double zcm[],double kx[],double ky[],double kz[],int seed,int N,double L,int n, int n_aglomerados,int n_aglomerados_x,int n_aglomerados_y,int n_aglomerados_z,double p_aglomerado,double theta_max, double Lx,double Ly, double Lz,double strain,double poisson);

void iteracion_mc(int mc,int first_seed,int N,double L,int n, int n_aglomerados,int n_aglomerados_x,int n_aglomerados_y,int n_aglomerados_z,double p_aglomerado,double theta_max,double d0,double dt, double Lx,double Ly, double Lz,double strain, double poisson,double R_pol,int percolado[],int numero_clusters_percolados[],int numero_cnt_percolando[], int numero_contactos_cnt_percolando[],double promedio_numero_contactos_cnt_percolando[],double d_prom_p[],double R[],int porcentaje_contactos);

//1.2	Main
int main(int argc, char *argv[])
{
	time_t ti,tf,ti_iter,tf_iter;
	ti = time(NULL);
		
	//se crea arreglo con los parametros que variarán en la simulación)
	int n_array[] = {1,5,9};
	double p_aglomerado_array[] = {0.0,0.1,0.2,0.3,0.4,0.6,0.8,1.0};
	double strain_array[] = {0.00,0.05,0.10,0.15,0.20};
	long long N_array[] = {1000,2000,3000,4000,5000,60000,8000,10000};
	
	int largo_n = sizeof(n_array)/sizeof(n_array[0]);
	int largo_p_aglomerado = sizeof(p_aglomerado_array)/sizeof(p_aglomerado_array[0]);
	int largo_strain = sizeof(strain_array)/sizeof(strain_array[0]);
	int largo_N = sizeof(N_array)/sizeof(N_array[0]);

	int nThreads = (argc > 1?atoi(argv[1]):0);
      int nIter,mc;
	std::thread t[nThreads];
		
	for (long long Ni=0; Ni<largo_N; Ni++){
		long long N = N_array[Ni];
	for (int n_i=0; n_i<largo_n; n_i++){
		int n = n_array[n_i];
	for (int Pi=0; Pi<largo_p_aglomerado; Pi++){
		double p_aglomerado = p_aglomerado_array[Pi];
	for (int strain_i=0; strain_i<largo_strain; strain_i++){
		double strain = strain_array[strain_i];

		ti_iter = time(NULL);
		double Lx_strain = Lx*(1+strain); //largo x del volumen representativo despues de la deformación
		double Ly_strain = Ly*(1-poisson*strain); //largo y del volumen representativo despues de la deformación
		double Lz_strain = Lz*(1-poisson*strain); //largo z del volumen representativo despues de la deformación
		double d_max_inicial_strain = sqrt(Lx_strain*Lx_strain + Ly_strain*Ly_strain + Lz_strain*Lz_strain);
		
		//se crean los arreglos en los que se guardaran los resultados para cada iteración de Monte Carlo
		int* percolado=new int[N_MonteCarlo](); //array que guarda la distancia promedio para cada iteracion
		int* numero_clusters_percolados=new int[N_MonteCarlo](); //array que guarda el numero de clusters percolados para cada iteracion
		int* numero_cnt_percolando=new int[N_MonteCarlo]();//array que guarda numero de cnt pertenecientes a un cluster percolado para cada iteracion
		int* numero_contactos_cnt_percolando=new int[N_MonteCarlo]();
		double* promedio_numero_contactos_cnt_percolando=new double[N_MonteCarlo]();
		double* d_prom_p=new double[N_MonteCarlo](); //array que guarda la distancia promedio para cada iteracion
		double* R=new double[N_MonteCarlo]();
	
		for (int ip=0; ip<N_MonteCarlo;++ip) //la distancia promedio se inicializa como la maxima distancia 
		{
			d_prom_p[ip] = d_max_inicial_strain;
			R[ip] = R_pol;
		}
		
	//se crean header del archivo si no existen
	if (guardar == 2){
		if (FILE *file = fopen("resultados_resistencia_v4_0/constantes.txt", "r")) {
	      fclose(file);
	    } else {
	      std::ofstream constantes_file;
			constantes_file.open("resultados_resistencia_v4_0/constantes.txt",
			std::ios::out | std::ios::app | std::ios::binary);     
			constantes_file<< "N L n nAglomerados nAglomeradosX nAglomeradosY nAglomeradosZ pAglomerado thetaMax d0 dvdw			<<" dt Lx Ly Lz strain poisson R_pol nMonteCarlo firstSeed"<<std::endl;    
			constantes_file.close();
	    } 
		
		if (FILE *file = fopen("resultados_resistencia_v4_0/resumen.txt", "r")) {
	      fclose(file);
	   } else {
	      std::ofstream resumen_file;
			resumen_file.open("resultados_resistencia_v4_0/resumen.txt",
std::ios::out | std::ios::app | std::ios::binary);
			resumen_file<<"N L n nAglomerados nAglomeradosX nAglomeradosY nAglomeradosZ pAglomerado thetaMax d0 dvdw dt
			<<" Lx Ly Lz strain poisson R_pol nMonteCarlo firstSeed percolado numeroCNTPercolando 
			<<" numeroContactosCNTPercolando PromedioNumeroContactosCNTPercolando R tiempo"<<std::endl; 
			resumen_file.close();  
	   } 
		
		if (FILE *file = fopen("resultados_resistencia_v4_0/resultados_por_seed.txt", "r")) {
	      fclose(file);
	   } else {
		std::ofstream resultados_por_seed_file;
		resultados_por_seed_file.open("resultados_resistencia_v4_0/resultados_por_seed.txt",
std::ios::out | std::ios::app | std::ios::binary);     
		resultados_por_seed_file<<"seed N L n nAglomerados nAglomeradosX nAglomeradosY nAglomeradosZ pAglomerado thetaMax d0 
		<<"dvdw dt Lx Ly Lz strain poisson R_pol percolado numeroCNTPercolando numeroContactosCNTPercolando 
		<<"PromedioNumeroContactosCNTPercolando porcentajeContactos numeroNodosA numeroElementosA promedioConexionesA R" 
		<<std::endl;  
		resultados_por_seed_file.close(); 
		}
	}	
		
	//se crean los arreglos en los que se guardaran los resultados despues de aplicar la deformacion para cada iteración de Monte Carlo
		mc = 0;
		do
		{
		int realT =0;
		for (nIter=0; nIter < nThreads && mc < N_MonteCarlo; nIter++)
			{
			// std::cout << "Lanzando thread " << nIter<<" en iteracion " << mc <<std::endl;

			t[nIter] = std::thread(iteracion_mc,mc,first_seed,N,L,n,n_aglomerados,n_aglomerados_x,n_aglomerados_y,n_aglomerados_z,p_aglomerado,theta_max,d0,dt,Lx,Ly,Lz,strain,poisson, R_pol,percolado,numero_clusters_percolados,numero_cnt_percolando,numero_contactos_cnt_percolando,promedio_numero_contactos_cnt_percolando,d_prom_p,R,porcentaje_contactos); //arreglos que se llenan en la funcion						
			mc++;
			realT++;
			}
		for (nIter=0; nIter < realT; nIter++)
			{
			// std::cout << "Esperando thread " << nIter <<std::endl;
			t[nIter].join();
			}
		} while (mc < N_MonteCarlo);
		
	tf_iter = time(NULL);
	double segundos_iter =   difftime(tf_iter,ti_iter);
	
	//se guardan los resultados en archivos .txt
	if (guardar == 2){	
		std::ofstream constantes_file;
		constantes_file.open("resultados_resistencia_v4_0/constantes.txt",
		std::ios::out | std::ios::app | std::ios::binary);     
		constantes_file<<std::fixed
		<<std::setprecision(0)<<N << " "<<L<< " "<<n<< " "<<n_aglomerados<< " "
		<<n_aglomerados_x<<" "<<n_aglomerados_y<<" "<<n_aglomerados_z<<" "
		<<std::setprecision(2)<<p_aglomerado<<" "
		<<std::setprecision(2) <<theta_max<<" "
		<<std::setprecision(0)<<d0<<" "<<dvdw<<" "<<dt<<" "<<Lx<<" "<<Ly<<" "<<Lz<<" "
		<<std::setprecision(2)<<strain<<" "<<poisson<<" "
		<<std::setprecision(0)<<R_pol<<" "
		<<std::setprecision(0)<<N_MonteCarlo<<" "<<first_seed<<" "
		<<std::endl;  
		constantes_file.close();
			
		std::ofstream numero_cnt_percolando_file;
		numero_cnt_percolando_file.open("resultados_resistencia_v4_0/numero_cnt_percolando.txt",
		std::ios::out | std::ios::app | std::ios::binary);     
		guardar1d(numero_cnt_percolando_file,numero_cnt_percolando,N_MonteCarlo);       
		numero_cnt_percolando_file.close();
		
		std::ofstream numero_contactos_cnt_percolando_file;
		numero_contactos_cnt_percolando_file.open("resultados_resistencia_v4_0/numero_contactos_cnt_percolando.txt",
		std::ios::out | std::ios::app | std::ios::binary);     
		guardar1d(numero_contactos_cnt_percolando_file,numero_contactos_cnt_percolando,N_MonteCarlo);       
		numero_contactos_cnt_percolando_file.close();		
		
		std::ofstream promedio_numero_contactos_cnt_percolando_file;
													 promedio_numero_contactos_cnt_percolando_file.open("resultados_resistencia_v4_0/
<<"promedio_numero_contactos_cnt_percolando.txt",std::ios::out | std::ios::app | std::ios::binary);     
		guardar1d(promedio_numero_contactos_cnt_percolando_file,promedio_numero_contactos_cnt_percolando,N_MonteCarlo);       
		promedio_numero_contactos_cnt_percolando_file.close();
		
		std::ofstream resumen_file;
		resumen_file.open("resultados_resistencia_v4_0/resumen.txt",
		std::ios::out | std::ios::app | std::ios::binary);     
		resumen_file<<std::fixed
		<<std::setprecision(0)<<N << " "<<L<< " "<<n<< " "<<n_aglomerados<< " "
		<<n_aglomerados_x<<" "<<n_aglomerados_y<<" "<<n_aglomerados_z<<" "
		<<std::setprecision(2)<<p_aglomerado<<" "
		<<std::setprecision(2) <<theta_max<<" "
		<<std::setprecision(0)<<d0<<" "<<dvdw<<" "<<dt<<" "<<Lx<<" "<<Ly<<" "<<Lz<<" "
		<<std::setprecision(2)<<strain<<" "<<poisson<<" "
		<<std::setprecision(0)<<R_pol<<" "
		<<std::setprecision(0)<<N_MonteCarlo<<" "<<first_seed<<" "
		<<std::setprecision(2)<<mean(percolado,N_MonteCarlo)<<" "
		<<std::setprecision(0)<<mean(numero_cnt_percolando,N_MonteCarlo)<<" "
		<<std::setprecision(0)<<mean(numero_contactos_cnt_percolando,N_MonteCarlo)<<" "
		<<std::setprecision(1)<<mean(promedio_numero_contactos_cnt_percolando,N_MonteCarlo)<<" "
		<<std::setprecision(5)<<mean(R,N_MonteCarlo)<<" "
		<<std::setprecision(0)<<segundos_iter
		<<std::endl;    
		resumen_file.close(); 
	}

	if (imprimir>=1)
	{
		std::cout <<std::fixed
		<<setprecision(0)<<"Box:"<< Lx/L 
		<<std::setprecision(0)<<" n:"<<n
		<<std::setprecision(0)<<" N:" <<N 
		<<std::setprecision(2)<<" p_a:" << p_aglomerado
		<<std::setprecision(0)<<" n_a:" << n_aglomerados
		<<std::setprecision(2)<<" str:" << strain 
		<<std::setprecision(0)<<" R_pol:" << R_pol 
		<<std::setprecision(2)<<" per:" <<mean(percolado,N_MonteCarlo)
		<<std::setprecision(0)<<"  n_cnt_perc:" <<mean(numero_cnt_percolando,N_MonteCarlo)
		<<std::setprecision(0)<<"  n_cont_cnt_perc:" <<mean(numero_contactos_cnt_percolando,N_MonteCarlo)
		<<std::setprecision(1)<<"  prom_cont_cnt_perc:" <<mean(promedio_numero_contactos_cnt_percolando,N_MonteCarlo)
		<<std::setprecision(5)<<" R:"<<mean(R,N_MonteCarlo)
		<<std::setprecision(0)<<"  t:"<<segundos_iter
		<<std::endl; 
		if (imprimir==2){
				std::cout <<"\n"<<std::endl;
		}
	}
	
	delete[] percolado;
	delete[] numero_clusters_percolados;
	delete[] numero_cnt_percolando;
	delete[] numero_contactos_cnt_percolando;
	delete[] promedio_numero_contactos_cnt_percolando;
	delete[] d_prom_p;
	delete[] R;
	
	} //end for N
	} //end for p_aglomerado
	} //end for strain
	} //end for n

 	
tf = time(NULL);	
   //se entrega el tiempo total utilizado en el cálculo de las iteraciones
double segundos_totales =   difftime(tf,ti);
std::cout << "tiempo total en segundos:" << segundos_totales << std::endl;  
return 0;
} //end main

//1.3	Funciones implementadas
/*funciones implementadas*/
int redondear(double x){
return std::ceil(x);}

double rand_xor128(int seed,long llamadas_por_seed) { //generador de numeros aleatorios en (0,1). En Matlab se implemetó el mismo generador
	
	long first_position = llamadas_por_seed*(seed-1)+1;
	static int last_seed = -1;
	static uint32_t x = 123456789;
	static uint32_t y = 362436069;
	static uint32_t z = 521288629;
	static uint32_t w = 88675123; 
	uint32_t t;
	
	if (last_seed != seed){

		if (first_position>=30000001){
			x = 577419830;
			y = 162715901;
			z = 1330228817;
			w = 501169124;
			for (int i=0;i<(first_position-1)-30000000;i++){
				t = x ^ (x << 11);
				x = y; y = z; z = w;
				w = w ^ (w >> 19) ^ (t ^ (t >> 8));		
			}			
		}	

		else if (first_position>=20000001){
			x = 3393161203;
			y = 4284861433;
			z = 2290513254;
			w = 2155665518;
			
			for (int i=0;i<(first_position-1)-20000000;i++){
				t = x ^ (x << 11);
				x = y; y = z; z = w;
				w = w ^ (w >> 19) ^ (t ^ (t >> 8));		
			}			
		}
		
		else if (first_position>=10000001){
			x = 2255198544;
			y = 1353851617;
			z = 1528002280;
			w = 2928528747;
			
			for (int i=0;i<(first_position-1)-10000000;i++){
				t = x ^ (x << 11);
				x = y; y = z; z = w;
				w = w ^ (w >> 19) ^ (t ^ (t >> 8));		
			}
		}
		
		else{ 
			x = 123456789;
			y = 362436069;
			z = 521288629;
			w = 88675123;
			
			for (int i=0;i<(first_position-1);i++){
				t = x ^ (x << 11);
				x = y; y = z; z = w;
				w = w ^ (w >> 19) ^ (t ^ (t >> 8));		
			}
		}
	}
	
	
	t = x ^ (x << 11);
	x = y; y = z; z = w;
	w = w ^ (w >> 19) ^ (t ^ (t >> 8));
	
	last_seed = seed;
	return (double)w/UINT32_MAX;
}

double dotp(double a[],double b[],int size){ //funcion que calcula el producto escalar entre dos array
	double suma = 0;
	for (int i=0;i<size;i++)
	{
		suma = suma + a[i]*b[i];
	}
	return suma;
}

double mean(int a[],int size){ //funcion que calcula el promedio de los elementos de un array tipo int
	double suma = 0;
	for (int i=0;i<size;i++)
	{
		suma = suma + a[i];
	}
	return (suma/size);
}	

double mean(double a[],int size){ //funcion que calcula el promedio de los elementos de un array tipo double
	double suma = 0;
	for (int i=0;i<size;i++)
	{
		suma = suma + a[i];
	}
	return (suma/size);
}

double sd(int a[],int size){ //funcion que calcula el promedio de los elementos de un array tipo int
	if (size == 1)
	{return -1;}
	else
	{
	double suma = 0;
	double avg = mean(a,size);
	for (int i=0;i<size;i++)
		{suma = suma + (a[i]-avg)*(a[i]-avg);}
	double var = suma/(size-1);
	return (sqrt(var));
	}
}	

double sd(double a[],int size){ //funcion que calcula el promedio de los elementos de un array tipo int
	double suma = 0;
	double avg = mean(a,size);
	for (int i=0;i<size;i++)
		{suma = suma + (a[i]-avg)*(a[i]-avg);	}
	double var = suma/(size-1);
	return (sqrt(var));
}

void imprimir1d(double myArray[], int size){ //funcion que imprime en pantalla un array de una dimension
	int i = 0;
	while (i < size)
	   {
	        cout<<setprecision(4)<<myArray[i]<<" "; //se imprime con 4 digitos de precision
	        ++i;
	   }
	cout << endl;
}

void imprimir1d(int myArray[], int size){ //funcion que imprime en pantalla un array de una dimension
	int i = 0;
	while (i < size)
	   {
	        cout<<setprecision(4)<<myArray[i]<<" "; //se imprime con 4 digitos de precision
	        ++i;
	   }
	cout << endl;
}

void imprimir2d(double matrix[], int rows, int cols){ //funcion que imprime en pantalla un array de dos dimensiones
   for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++)
            {cout << setprecision(5)<<matrix[i*cols+j]<<" ";}
        cout << endl;
   }
}

void imprimir2d(int matrix[], int rows, int cols){ //funcion que imprime en pantalla un array de dos dimensiones
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++)
            {cout << setprecision(5)<<matrix[i*cols+j]<<" ";}
        cout << endl;
    }
}

void imprimir2dmin(double matrix[], int rows, int cols){ //funcion que imprime en pantalla un array de dos dimensiones
   for (int i = 0; i < rows; i++){
		cout << setprecision(5)<<matrix[i*cols+1]+1<<" ";
		cout << setprecision(5)<<matrix[i*cols+0]+1<<" ";
		cout << setprecision(5)<<matrix[i*cols+2]<<" ";
        cout << endl;
   }
}

void guardar1d(std::ofstream &outputFile,int myArray[], int size){ //funcion que guarda en archivo un array tipo int de una dimension
	for (int i = 0; i < size; i++)
        {outputFile << myArray[i]<<" ";}
    outputFile << endl;
}

void guardar1d(std::ofstream &outputFile,double myArray[], int size){ //funcion que guarda en archivo un array tipo double de una dimension
for (int i = 0; i < size; i++)
        {outputFile << myArray[i]<<" ";}
    outputFile << endl;
}

void guardar2d(std::ofstream &outputFile,double matrix[], int rows, int cols){ //funcion que guarda en archivo un array de dos dimensiones
    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++)
            {outputFile << matrix[i*cols+j]<<" ";}
        outputFile << endl;
    }
}

void guardar_constantes(std::ofstream &outputFile,int seed,int N,double L,int n, int n_aglomerados,int n_aglomerados_x,int n_aglomerados_y,int n_aglomerados_z,
				double p_aglomerado,double theta_max,double d0,double dt, double Lx,double Ly, double Lz,double strain, double poisson){ 
				//funcion que guarda en archivo un array de dos dimensiones
	outputFile <<std::fixed
	<<std::setprecision(0)<<seed<<" "<<N <<" "<<L<<" "<<n<<" "<<n_aglomerados<<" "
	<<n_aglomerados_x<<" "<<n_aglomerados_y<<" "<<n_aglomerados_z<<" "
	<<std::setprecision(2)<<p_aglomerado<<" "
	<<std::setprecision(2) <<theta_max<<" "
	<<std::setprecision(0)<<d0<<" "<<dvdw<<" "<<dt<<" "<<Lx<<" "<<Ly<<" "<<Lz<<" "
	<<std::setprecision(2)<<strain<<" "<<poisson<<" "
	<<std::endl;  
} 

double G_tunel(double d,double lambda){ //funcion que calcula la resistencia por efecto tunel a partir de la distancia entre dos cnt
	// double h = 6.6261e-34; //[Js]
	// double e = 1.6022e-19; //[C]
	// double m = 9.109e-31; //[kg]
	// double A = (d0*1e-9)*(d0*1e-9); //[m^2]
	// double lambda = 0.2 * e; //[J]
	double alpha = 31588.2957554212; //alpha = e*e*sqrt(2*m*e)/(h*h);
	double beta = 10.2461471970158; //beta = 4*pi*1e-9*sqrt(2*m*e)/h;
	double g = alpha*sqrt(lambda)/(exp(beta*d*sqrt(lambda)) );
	
	return g;
}


void bubble_sortA(double A_suma[], int numero_nodos) {
      bool swapped = true;
      int j = 0;
      int tmpi;
	  int tmpj;
	  double tmpx;
      while (swapped) {
            swapped = false;
            j++;
            for (int i = 0; i < numero_nodos - j; i++) {
                  if (A_suma[i*3+0] > A_suma[(i + 1)*3+0]) {
                     tmpj = A_suma[i*3+1];
                     tmpi = A_suma[i*3+0];
							tmpx = A_suma[i*3+2];
                     A_suma[i*3+1] = A_suma[(i + 1)*3+1];
                  	A_suma[i*3+0] = A_suma[(i + 1)*3+0];
                  	A_suma[i*3+2] = A_suma[(i + 1)*3+2];
                     A_suma[(i + 1)*3+1] = tmpj;
                     A_suma[(i + 1)*3+0] = tmpi;
                     A_suma[(i + 1)*3+2] = tmpx;
                     swapped = true;
                  }
            }
      }
}
	
void spice(std::ofstream &outputFile,double dmin_array[], int length_bb,int nodo_extremo[],int N,double lambda){ //Se crea archivo para ser ingresado a SPICE
	outputFile <<"*Archivo Spice resistencias tunel" << endl;   
    for (int k=0;k<length_bb;k++)
	{		
		int i = dmin_array[k*3];
		int j = dmin_array[k*3+1];
		double R = 1/G_tunel(dmin_array[k*3+2],lambda);
		if (nodo_extremo[j]==-1)
		{outputFile << "R" <<i+1<<"o"<<j+1<<" "<<i+1<<" "<<0<<" "<<R<<endl;}
		else if (nodo_extremo[j]==1)
		{outputFile << "R" <<i+1<<"o"<<j+1<<" "<<i+1<<" "<<N+1<<" "<<R<<endl;}
		else if (nodo_extremo[i]==-1)
		{outputFile << "R" <<i+1<<"o"<<j+1<<" "<<0<<" "<<j+1<<" "<<R<<endl;}
		else if (nodo_extremo[i]==1)
		{outputFile << "R" <<i+1<<"o"<<j+1<<" "<<N+1<<" "<<j+1<<" "<<R<<endl;}
		else
		{outputFile << "R" <<i+1<<"o"<<j+1<<" "<<i+1<<" "<<j+1<<" "<<R<<endl;}
	} 

outputFile << "V1 "<< 0 <<" "<< N+1 <<" "<< 1 << endl;
outputFile << ".op" << endl;
outputFile << ".PRINT I(V1)" << endl;
outputFile << ".END" << endl;
}

//1.3.1	Función configuración inicial
void configuracion_inicial(double x[], double y[],double z[],double theta_suma[],double phi_suma[],double xcm[],double ycm[],double zcm[],double kx[],double ky[],double kz[],int seed,int N,double L,int n, int n_aglomerados,int n_aglomerados_x,int n_aglomerados_y,int n_aglomerados_z,double p_aglomerado,double theta_max, double Lx,double Ly, double Lz,double strain,double poisson)
	{
	std::lock_guard<std::mutex> guard(myMutex);
	double l = L/n;
	double Lx_strain = Lx*(1+strain); //largo x del volumen representativo despues de la deformación
	double Ly_strain = Ly*(1-poisson*strain); //largo y del volumen representativo despues de la deformación
	double Lz_strain = Lz*(1-poisson*strain); //largo z del volumen representativo despues de la deformación
	double* coef_aglomerado=new double[N]();
	int* aglomerado=new int [N]() ;
	double* Lxc=new double[N]();	
	double* Lyc=new double[N]();	
	double* Lzc=new double[N]();
	double* x_aglomerados_fijos=new double[n_aglomerados_x]();
	double* y_aglomerados_fijos=new double[n_aglomerados_y]();
	double* z_aglomerados_fijos=new double[n_aglomerados_z]();
	double* x_aglomerados_aleatorios=new double[n_aglomerados]();
	double* y_aglomerados_aleatorios=new double[n_aglomerados]();
	double* z_aglomerados_aleatorios=new double[n_aglomerados]();
	double* aglomerado_al_que_pertence_x=new double[N]();
	double* aglomerado_al_que_pertence_y=new double[N]();
	double* aglomerado_al_que_pertence_z=new double[N]();
	double R_sigma;
	double sigma;
	double* x_cnt_aglomerado=new double[N]();	
	double* y_cnt_aglomerado=new double[N]();	
	double* z_cnt_aglomerado=new double[N]();
	double* x_cnt_no_aglomerado=new double[N]();	
	double* y_cnt_no_aglomerado=new double[N]();	
	double* z_cnt_no_aglomerado=new double[N]();
	double* Ux1=new double[N]();	
	double* Uy1=new double[N]();	
	double* Uz1=new double[N]();
	double* Ux2=new double[N]();	
	double* Uy2=new double[N]();	
	double* Uz2=new double[N]();		
	
	long llamadas_por_seed = N*(13+2*n)+3*n_aglomerados;

	for (int i=0;i<N;i++){ //se establece si los cnt estan aglomerados o no de acuerdo a p_aglomerado
	coef_aglomerado[i] = rand_xor128(seed,llamadas_por_seed);
		if (coef_aglomerado[i]<p_aglomerado)
			{aglomerado[i]=1;}
		else
			{aglomerado[i]=0;}	
	}
	
	/*se generan la posiciones de los aglomerados y se asignan los cnt a cada uno de ellos*/
	if (n_aglomerados == 0){ //caso numero fijo de aglomerados distribuidos uniformemente en cada direccion 
		for (int i=0;i<n_aglomerados_x;i++) //posicion en x de los centros de los aglomerados 
		{x_aglomerados_fijos[i] = Lx*((double)(2*(i+1)-1)/(2*n_aglomerados_x));}
	
		for (int i=0;i<n_aglomerados_y;i++) //posicion en y de los centros de los aglomerados 
		{y_aglomerados_fijos[i] = Ly*((double)(2*(i+1)-1)/(2*n_aglomerados_y));}
		
		for (int i=0;i<n_aglomerados_z;i++) //posicion en z de los centros de los aglomerados 
		{z_aglomerados_fijos[i] = Lz*((double)(2*(i+1)-1)/(2*n_aglomerados_z));}
		
		// aglomerado al que pertenece cada uno de los N nanotubos. Se hace un "for" separados para generar numeros igual a Matlab
		for (int i=0;i<N;i++)
			{aglomerado_al_que_pertence_x[i] = redondear(n_aglomerados_x*rand_xor128(seed,llamadas_por_seed));}
		for (int i=0;i<N;i++)
			{aglomerado_al_que_pertence_y[i] = redondear(n_aglomerados_y*rand_xor128(seed,llamadas_por_seed));}
		for (int i=0;i<N;i++)
			{aglomerado_al_que_pertence_z[i] = redondear(n_aglomerados_z*rand_xor128(seed,llamadas_por_seed));}	    
		
	    for(int i=0;i<N;++i) { //cada cnt se asocia a la posicion central de uno de los aglomerados
	        int j = aglomerado_al_que_pertence_x[i]-1;
			Lxc[i] = x_aglomerados_fijos[j];
			j = aglomerado_al_que_pertence_y[i]-1;
	        Lyc[i] = y_aglomerados_fijos[j];    
	        j = aglomerado_al_que_pertence_z[i]-1;
	        Lzc[i] = z_aglomerados_fijos[j];
	    }
	} //end if n_aglomerado==0      
	
	else{  //caso aglomerados posicionados aleatoriamente
		//posicion en x,y,z de los centros de los aglomerados, generados aleatoriamente
		for (int i=0;i<n_aglomerados;i++)
		{x_aglomerados_aleatorios[i] = Lx*rand_xor128(seed,llamadas_por_seed);}
		for (int i=0;i<n_aglomerados;i++)
		{y_aglomerados_aleatorios[i] = Ly*rand_xor128(seed,llamadas_por_seed);}
		for (int i=0;i<n_aglomerados;i++)
		{z_aglomerados_aleatorios[i] = Lz*rand_xor128(seed,llamadas_por_seed);}
		
		// aglomerado al que pertenece cada uno de los N nanotubos. Se hace un "for" separados para generar numeros igual a Matlab
		for (int i=0;i<N;i++)
			{aglomerado_al_que_pertence_x[i] = redondear(n_aglomerados*rand_xor128(seed,llamadas_por_seed));}
		for (int i=0;i<N;i++)
			{aglomerado_al_que_pertence_y[i] = redondear(n_aglomerados*rand_xor128(seed,llamadas_por_seed));}
		for (int i=0;i<N;i++)
			{aglomerado_al_que_pertence_z[i] = redondear(n_aglomerados*rand_xor128(seed,llamadas_por_seed));}
	
		for(int i=0;i<N;++i){ //cada cnt se asocia a laposicion central de uno de los aglomerados
	        int j = aglomerado_al_que_pertence_x[i]-1;
			Lxc[i] = x_aglomerados_aleatorios[j];
			j = aglomerado_al_que_pertence_y[i]-1;
	        Lyc[i] = y_aglomerados_aleatorios[j];    
	        j = aglomerado_al_que_pertence_z[i]-1;
	        Lzc[i] = z_aglomerados_aleatorios[j];
	    }
	} //end else
	
/*se generan las posiciones de los nanotubos aglomerados siguiento una distribución normal en torno a
los centros de los aglomerados. La distribución normal se logra usando el metodo de Box-Muller*/
	//se calcula la desviación estandar de modo que el radio este dentro del volumen respresentativo
	if (n_aglomerados == 0)
		{R_sigma = 1/(double)(2*std::max(n_aglomerados_x,std::max(n_aglomerados_y,n_aglomerados_z)));}
	else
		{R_sigma = (double)1/(2*std::cbrt((double)n_aglomerados));}

	sigma = R_sigma/sqrt(2*log(N)); //R = sqrt(2*log(n))*SD
	
	//se generan numeros aleatorios para generar la distribución segun Box-Muller
	for (int i=0;i<N;i++)		
	{Ux1[i] = rand_xor128(seed,llamadas_por_seed);}
	for (int i=0;i<N;i++)		
	{Ux2[i] = rand_xor128(seed,llamadas_por_seed);}		
	for (int i=0;i<N;i++)		
	{Uy1[i] = rand_xor128(seed,llamadas_por_seed);}
	for (int i=0;i<N;i++)		
	{Uy2[i] = rand_xor128(seed,llamadas_por_seed);}
	for (int i=0;i<N;i++)		
	{Uz1[i] = rand_xor128(seed,llamadas_por_seed);}
	for (int i=0;i<N;i++)		
	{Uz2[i] = rand_xor128(seed,llamadas_por_seed);}
	
	//Metodo Box-Muller                           
	for (int i=0;i<N;i++){
		x_cnt_aglomerado[i] = Lxc[i] + Lx*sigma*sqrt(-2*log(Ux1[i]))*cos(2*pi*Ux2[i]);
		y_cnt_aglomerado[i] = Lyc[i] + Ly*sigma*sqrt(-2*log(Uy1[i]))*cos(2*pi*Uy2[i]);
		z_cnt_aglomerado[i] = Lzc[i] + Lz*sigma*sqrt(-2*log(Uz1[i]))*cos(2*pi*Uz2[i]);
	}
	
	//se generan las posiciones de los nanotubos no aglomerados al azar
	for (int i=0;i<N;i++)
		{x_cnt_no_aglomerado[i] = Lx*rand_xor128(seed,llamadas_por_seed);}
	for (int i=0;i<N;i++)
		{y_cnt_no_aglomerado[i] = Ly*rand_xor128(seed,llamadas_por_seed);}
	for (int i=0;i<N;i++)
		{z_cnt_no_aglomerado[i] = Lz*rand_xor128(seed,llamadas_por_seed);}
	
	//Se calcula la posicion del centro de masa de cada cnt. Se considera si esta aglomerado o no
	int ncm = (n-1)/2;
	for (int i=0;i<N;i++){
		xcm[i*n +ncm] = x_cnt_aglomerado[i]*aglomerado[i] + x_cnt_no_aglomerado[i]*(1-aglomerado[i]);
		ycm[i*n +ncm] = y_cnt_aglomerado[i]*aglomerado[i] + y_cnt_no_aglomerado[i]*(1-aglomerado[i]);
		zcm[i*n +ncm] = z_cnt_aglomerado[i]*aglomerado[i] + z_cnt_no_aglomerado[i]*(1-aglomerado[i]);
	}
	
	//angulos entre los segmentos de cada nanotubo
	for (int j=0;j<n;j++){ //theta en (-theta_max,theta_max)
		for (int i=0;i<N;i++)
			{theta_suma[i*n+j] = acos((theta_max/(pi/2))*(rand_xor128(seed,llamadas_por_seed)*2-1));}//theta en [pi/2 -theta_max, pi/2 + theta_max]
	}
	for (int j=0;j<n;j++){ //phi en (-theta_max,theta_max)
		for (int i=0;i<N;i++)
		{phi_suma[i*n+j] = (rand_xor128(seed,llamadas_por_seed)*2-1)*theta_max;} //phi en [-theta_max, theta_max]
	}
	
	//DEFORMACION
	if (strain>0){
		for (int i=0;i<N;i++){
			xcm[i*n +ncm] = xcm[i*n +ncm]*(1+strain);
			ycm[i*n +ncm] = ycm[i*n +ncm]*(1-poisson*strain);
			zcm[i*n +ncm] = zcm[i*n +ncm]*(1-poisson*strain);
		}
	}
	
	double y_atan,x_atan;
	double* theta_suma_s = new double [N*n];
	double* phi_suma_s = new double [N*n];
	if (strain>0){
		for (int j=0;j<n;j++) //phi en (-theta_max,theta_max)
		{
			for (int i=0;i<N;i++)
			{
				y_atan = sin(phi_suma[i*n+j])*(1-poisson*strain)/(1+strain);
				x_atan = cos(phi_suma[i*n+j]);
				phi_suma_s[i*n+j] = atan2(y_atan,x_atan);
			}
		}	
		
		for (int j=0;j<n;j++) //theta en (-theta_max,theta_max)
		{
			for (int i=0;i<N;i++)
				{
					y_atan = sin(theta_suma[i*n+j])*sin(phi_suma[i*n+j])/sin(phi_suma_s[i*n+j]);
					x_atan = cos(theta_suma[i*n+j]);
					theta_suma_s[i*n+j] = atan2(y_atan,x_atan);
				}
		}
		
		for (int j=0;j<n;j++) //reasignar los nuevos valores a phi_suma y theta_suma
		{
			for (int i=0;i<N;i++)
			{
				phi_suma[i*n+j] = phi_suma_s[i*n+j]; 
				theta_suma[i*n+j] = fmod(theta_suma_s[i*n+j]+2*pi,2*pi); 
			}
		}	
	}
	
	//Se calculan las posiciones de las extremo del segmento central del cnt
	for(int i=0;i<N;i++){ 
		//extremo izquierdo
		x[i*(n+1)+ncm] = xcm[i*n+ncm] - (l/2)*sin(theta_suma[i*n+ncm])*cos(phi_suma[i*n+ncm]); 
		y[i*(n+1)+ncm] = ycm[i*n+ncm] - (l/2)*sin(theta_suma[i*n+ncm])*sin(phi_suma[i*n+ncm]);
		z[i*(n+1)+ncm] = zcm[i*n+ncm] - (l/2)*cos(theta_suma[i*n+ncm]);
		
		//extremo derecho
		x[i*(n+1)+(ncm+1)] = xcm[i*n+ncm] + (l/2)*sin(theta_suma[i*n+ncm])*cos(phi_suma[i*n+ncm]); 
		y[i*(n+1)+(ncm+1)] = ycm[i*n+ncm] + (l/2)*sin(theta_suma[i*n+ncm])*sin(phi_suma[i*n+ncm]);
		z[i*(n+1)+(ncm+1)] = zcm[i*n+ncm] + (l/2)*cos(theta_suma[i*n+ncm]);
	}

	//Se calculan las posiciones de los siguientes extremos
	for(int i=0;i<N;i++){ 
		//extremos hacia la izquierda
		for(int j=ncm;j>0;j--){
			x[i*(n+1)+j-1] = x[i*(n+1)+j] - l*sin(theta_suma[i*n+(j-1)])*cos(phi_suma[i*n+(j-1)]); 
			y[i*(n+1)+j-1] = y[i*(n+1)+j] - l*sin(theta_suma[i*n+(j-1)])*sin(phi_suma[i*n+(j-1)]);
			z[i*(n+1)+j-1] = z[i*(n+1)+j] - l*cos(theta_suma[i*n+(j-1)]);
		}
		
		//extremos hacia la derecha
		for(int j=ncm+1;j<n;j++){
			x[i*(n+1)+j+1] = x[i*(n+1)+j] + l*sin(theta_suma[i*n+j])*cos(phi_suma[i*n+j]); 
			y[i*(n+1)+j+1] = y[i*(n+1)+j] + l*sin(theta_suma[i*n+j])*sin(phi_suma[i*n+j]);
			z[i*(n+1)+j+1] = z[i*(n+1)+j] + l*cos(theta_suma[i*n+j]);
		}
	}
	
//	cout<<"x:\n"; imprimir2d((double*)x,N,n+1); cout<<endl;
//	cout<<"y:\n"; imprimir2d((double*)y,N,n+1); cout<<endl;
//	cout<<"z:\n"; imprimir2d((double*)z,N,n+1); cout<<endl;
	
	//se calculan las posiciones de los centros de masa de cada segmento
	for(int i=0;i<N;i++){ 
		for(int j=0;j<n;j++){
			xcm[i*n+j] = (x[i*(n+1)+j] + x[i*(n+1)+j+1])/2;
			ycm[i*n+j] = (y[i*(n+1)+j] + y[i*(n+1)+j+1])/2;
			zcm[i*n+j] = (z[i*(n+1)+j] + z[i*(n+1)+j+1])/2;
		}
	}
	
	//se calculan los vectores unitarios asociados a la orientacion de cada segmento
	for(int i=0;i<N;i++){ 
		for(int j=0;j<n;j++){
			kx[i*n+j] = sin(theta_suma[i*n+j])*cos(phi_suma[i*n+j]); 
			ky[i*n+j] = sin(theta_suma[i*n+j])*sin(phi_suma[i*n+j]);
			kz[i*n+j] = cos(theta_suma[i*n+j]);
		}
	}
	
//	cout<<"kx:\n"; imprimir2d((double*)kx,N,n); cout<<endl;
//	cout<<"ky:\n"; imprimir2d((double*)ky,N,n); cout<<endl;
// cout<<"kz:\n"; imprimir2d((double*)kz,N,n); cout<<endl;

	delete[] coef_aglomerado;
	delete[] aglomerado;
	delete[] Lxc;	
	delete[] Lyc;	
	delete[] Lzc;
	delete[] x_aglomerados_fijos;
	delete[] y_aglomerados_fijos;
	delete[] z_aglomerados_fijos;
	delete[] x_aglomerados_aleatorios;
	delete[] y_aglomerados_aleatorios;
	delete[] z_aglomerados_aleatorios;
	delete[] aglomerado_al_que_pertence_x;	
	delete[] aglomerado_al_que_pertence_y;	
	delete[] aglomerado_al_que_pertence_z;
	delete[] x_cnt_aglomerado;	
	delete[] y_cnt_aglomerado;	
	delete[] z_cnt_aglomerado;
	delete[] x_cnt_no_aglomerado;	
	delete[] y_cnt_no_aglomerado;	
	delete[] z_cnt_no_aglomerado;
	delete[] Ux1;	
	delete[] Ux2;	
	delete[] Uy1;	
	delete[] Uy2;	
	delete[] Uz1;	
	delete[] Uz2;
	
	rand_xor128(-2,-1); //se llama rand con un seed negativo para reiniciar el generador
} //end configuracion_inicial

//1.3.2	Función iteración de Monte Carlo
void iteracion_mc(int mc,int first_seed,int N,double L,int n, int n_aglomerados,int n_aglomerados_x,int n_aglomerados_y,int n_aglomerados_z,double p_aglomerado,double theta_max,double d0,double dt, double Lx,double Ly, double Lz,double strain,double poisson, double R_pol,int percolado[],int numero_clusters_percolados[],int numero_cnt_percolando[], int numero_contactos_cnt_percolando[],double promedio_numero_contactos_cnt_percolando[],double d_prom_p[],double R[],int porcentaje_contactos)
{	
	double l = L/n;
	double Lx_strain = Lx*(1+strain); //largo x del volumen representativo despues de la deformación
	double Ly_strain = Ly*(1-poisson*strain); //largo y del volumen representativo despues de la deformación
	double Lz_strain = Lz*(1-poisson*strain); //largo z del volumen representativo despues de la deformación
	double d_max_inicial_strain = sqrt(Lx_strain*Lx_strain + Ly_strain*Ly_strain + Lz_strain*Lz_strain);
	int seed = first_seed+mc;
			  	
	double* x=new double[N*(n+1)]();	
	double* y=new double[N*(n+1)]();	
	double* z=new double[N*(n+1)]();
	double* theta_suma=new double[N*n]();	
	double* phi_suma=new double[N*n]();	
	double* xcm=new double[N*n]();	
	double* ycm=new double[N*n]();	
	double* zcm=new double[N*n]();
	double* kx=new double[N*n]();	
	double* ky=new double[N*n]();	
	double* kz=new double[N*n]();
	double dcm = d_max_inicial_strain;
	int n_2 = (n-1)/2;
	double xcm_2 = Lx_strain;
	double ycm_2 = Ly_strain;
	double zcm_2 = Lz_strain;
	double* rcm_i=new double[3*n]();	
	double* rcm_j=new double[3*n]();	
	double* dcmkm=new double[n*n]();
	double* rki=new double[3]();	
	double* rmj=new double[3]();	
	double dcmkm_min = d_max_inicial_strain;
	double* ri=new double[3]();	
	double* rj=new double[3]();	
	double* p_i=new double[3]();
	double* pj=new double[3]();
	double* rij=new double[3]();
	double* rji=new double[3]();
	double* v = new double[3]();	
	double dotp_pi_pj = 0;
	double sij = 0;
	double sji = 0;
	double l_2 = l/2;
	int num_cluster = -1;
	double xmin = Lx_strain;
	double xmax = 0;
	int* cluster=new int[N]();
	for (int i=0;i<N;i++)
	{cluster[i] = i+1;}
	int* indice_cluster_percolado_aux=new int[N]();
	int node_min_x = -1;
	int node_max_x = -1;
	int node_min_lista = -1;
	int node_max_lista = -1;
	int node_min_tunel = -1;
	int node_max_tunel = -1;
	int numero_contactos = 0;
	int numero_posibles_contactos = 0;
	int numero_posibles_contactos_segmentos = 0;
	int numero_contactos_penetrados = 0;	
	int* cnt_percolado=new int[N]();
	long long N_5 = N/5;
	long long largo_dmin = (N_5*N_5/4*3*porcentaje_contactos);
	double* dmin=new double[largo_dmin]();	
	int* numero_contactos_cnt=new int[N]();	
	int* nodo_extremo = new int[N]();
	double* dmin_bb=new double[largo_dmin]();
	int* nodo_bb = new int[N]();	
	long numero_edges_bb = 0;
	double dmin_aux = d_max_inicial_strain;
	double* A_pre=new double[4*largo_dmin]();	
	int* nodo_de_A = new int[N+1]();
	int* indice_cnt_A = new int[2*N]();
	// int* indice_cnt_dmin_bb = new int[2*N]();
	
	int cambiados = 1;
	long indice_dmin = 0;
	int eliminados = 1;
	int primer_nodo_max = -1;
	int primer_nodo_min = 1;
	int ipnm = 0;
	int indice_A_pre = 0;	
	int numero_elementos_A = 0;
	int numero_nodos_A = 0;
	int contador_Ap = 0;
	double V = 1;
	double I = -0.0001;
			
	for(int k=0;k<largo_dmin/3;k++){
				dmin[k] = -1;
				dmin[k+1] = -1;
				dmin[k+2] = d_max_inicial_strain;
		}	

	configuracion_inicial(x,y,z,theta_suma,phi_suma,xcm,ycm,zcm,kx,ky,kz,seed,N,L,n,n_aglomerados,n_aglomerados_x,
												n_aglomerados_y,n_aglomerados_z,p_aglomerado,theta_max,Lx,Ly,Lz,strain,poisson);	

	/* Cálculo de distancia minima entre los CNT*/		
	for (int i = 0;i<N;++i){ //calculo de la distancia mínima entre los nanotubos i,j
	    for (int j=i+1;j<N;++j){ 
	    dmin_aux = d_max_inicial_strain;
			xcm_2 = (xcm[i*n+n_2]-xcm[j*n+n_2]);
			ycm_2 = (ycm[i*n+n_2]-ycm[j*n+n_2]);
			zcm_2 = (zcm[i*n+n_2]-zcm[j*n+n_2]);
			
			dcm = sqrt(xcm_2*xcm_2 + ycm_2*ycm_2 + zcm_2*zcm_2);
			if (dcm < L + d0 + dt){
				for (int s=0;s<n;s++){
					rcm_i[0*n+s] = xcm[i*n+s];
					rcm_i[1*n+s] = ycm[i*n+s];
					rcm_i[2*n+s] = zcm[i*n+s];
				
					rcm_j[0*n+s] = xcm[j*n+s];
					rcm_j[1*n+s] = ycm[j*n+s];
					rcm_j[2*n+s] = zcm[j*n+s];
				}
					
				for(int k=0; k<n; ++k){//distancia del cm del fragmento k del nanotubo i al cm del fragmento m del nanotubo j
					for(int m=0; m<n; ++m){  
						for(int s=0;s<3;s++){
						rki[s] = rcm_i[s*n+k];
						rmj[s] = rcm_j[s*n+m];       
						}		
						dcmkm[k*n+m] = sqrt((rki[0]-rmj[0])*(rki[0]-rmj[0])+
							(rki[1]-rmj[1])*(rki[1]-rmj[1])+(rki[2]-rmj[2])*(rki[2]-rmj[2]));
					}//end for m
				} //end for k
			
			/*si la distancia entre los cm de los segmentos permite efecto tunel se calcula la distancia mínima entre dichos segmentos*/
				
				double dmin_ab = d_max_inicial_strain;				    
				
				for (int a=0; a<n; a++){
					for (int b=0; b<n; b++){
						if (  ( (dcmkm[a*n+b]<(l+d0+1*dt)) && (dvdw>0) ) 
							 || ( (dcmkm[a*n+b]<(l+1*dt)) && (dvdw == 0) )  ){ 
							for(int s=0;s<3;s++)
							{
								ri[s] = rcm_i[s*n+a];
								rj[s] = rcm_j[s*n+b];
							}
								p_i[0] = kx[i*n+a]; p_i[1] = ky[i*n+a]; p_i[2] = kz[i*n+a];
								pj[0] = kx[j*n+b]; pj[1] = ky[j*n+b]; pj[2] = kz[j*n+b];
							for(int s=0;s<3;s++)
							{
								rij[s] = ri[s]-rj[s];
								rji[s] = rj[s]-ri[s];
							}
							dotp_pi_pj = dotp(p_i,pj,3);
							
							if ( (p_i[0] != pj[0]) && (p_i[1] != pj[1]) && (p_i[2] != pj[2]))
							{
								sij = (dotp(rij,pj,3)*dotp_pi_pj-dotp(rij,p_i,3))
											/(1-(dotp_pi_pj)*(dotp_pi_pj));         
								sji = (dotp(rji,p_i,3)*dotp_pi_pj-dotp(rji,pj,3))
											/(1-(dotp_pi_pj)*(dotp_pi_pj));
							}
							
							if (abs(sij)<l_2 && abs(sji)<l_2){ //interaccion lado-lado
								sij = sij;
								sji = sji;
							}
							else if (abs(sij)<l_2 && sji>l_2){ //interaccion extremo-lado
								for (int s=0;s<3;s++)
								{v[s] = rj[s] + l_2*pj[s] - ri[s];}
								sij = dotp(p_i,v,3);
								sji = l_2;
							}
							else if (abs(sij)<l_2 && sji<-l_2){ //interaccion extremo-lado
								for (int s=0;s<3;s++)
								{v[s] = rj[s] - l_2*pj[s] - ri[s];}
								sij = dotp(p_i,v,3);
								sji = -l_2;
							}
							else if (sij>l_2 && abs(sji)<l_2){ //interaccion extremo-lado
								for (int s=0;s<3;s++)
								{v[s] = ri[s] + l_2*p_i[s] - rj[s];}
								sij = l_2;
								sji = dotp(pj,v,3);
							}
							else if ((sij < -l_2) && abs(sji)<l_2){ //interaccion extremo-lado
								for (int s=0;s<3;s++)
								{v[s] = ri[s] - l_2*p_i[s] - rj[s];}
								sij = -l_2;
								sji = dotp(pj,v,3);
							}
							else if (sij>l_2 && sji>l_2){ //interaccion extremo-extremo
								sij = l_2;
								sji = l_2;
							}
							else if (sij>l_2 && sji<-l_2){ //interaccion extremo-extremo
								sij = l_2;
								sji = -l_2;
							}
							else if (sij<-l_2 && sji>l_2){ //interaccion extremo-extremo
								sij = -l_2;
								sji = l_2;
							}
							else if (sij<-l_2 && sji<-l_2){ //interaccion extremo-extremo
								sij = -l_2;
								sji = -l_2;
							}
							else{ //nanotubos paralelos
								sij = 0;
								sji = 0;
							}

							for (int s=0;s<3;s++)
							{v[s] = rj[s] + sji*pj[s] - ri[s] - sij*p_i[s];}
							dmin_ab = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
							// dmin_ab = dmin_ab - d0;
							
							if(dmin_ab<=dmin_aux){
								if (dmin_aux == d_max_inicial_strain) //solo se gurda posible contacto entre nanotubos si es la primera vez que escuentra distancia menor
									{numero_posibles_contactos +=1;}
								dmin_aux = dmin_ab;
							}
							numero_posibles_contactos_segmentos += 1;
						} //end if dcmkm posible contacto	
					} //end for b
				} //end for a	
				
				if (dvdw > 0){
					dmin_aux = dmin_aux - d0;
					if (dmin_aux<dvdw){
						dmin_aux = dvdw;
						numero_contactos_penetrados += 1;
					}
				}
				
				/*fin del calculo de la distancia mínima entre los segmentos mas cercanos*/
				
				if (dmin_aux < dt){ //si existe efecto tunel entre los nanotubos, se les asigna el mismo cluster
					num_cluster = std::min(cluster[i],cluster[j]);
					cluster[i] = num_cluster;
					cluster[j] = num_cluster;
					dmin[numero_contactos*3] = i;
					dmin[numero_contactos*3+1] = j;
					dmin[numero_contactos*3+2] = dmin_aux;
					numero_contactos += 1;
					if (numero_contactos > (largo_dmin/3)){
						cout<< "Se supero capacidad de dmin. El numero de contactos es "
						<<numero_contactos<< ". Largo dmin es "<< largo_dmin<<endl;
						return;
					}						
				} // end if tunel
			} // end if dcm < L +d0 +dt	
	    }//end for j	
	} //end for i
	// cout<<"\n dmin:"<<endl;	imprimir2dmin(dmin,numero_contactos,3);
	
	// SE REPASA LA ASIGNACION DE CLUSTERS
	cambiados = 1;
	while (cambiados>0){
		cambiados = 0;
		for (int k=0;k<numero_contactos;k++){
			int i = dmin[k*3];
			int j = dmin[k*3+1];
			num_cluster = std::min(cluster[i],cluster[j]);
			if ( cluster[i]>num_cluster || cluster[j]>num_cluster )
				{
				cluster[i] = num_cluster;
				cluster[j] = num_cluster;
				cambiados +=1;
				}	
		} // end for k
	} // end while
	// cout<<"cluster:\n"<<endl;	imprimir1d(cluster,N);

	/*CALCULO DE RESULTADOS Y VARIABLES REPRESENTATIVAS*/
	// Determinar si existe un cluster que percole
	for (int i=0;i<N;++i){ 
		xmin = Lx_strain;
		xmax = 0;
	    for (int j=0;j<N;j++){ // calculo de las posiciones mínimas y máximas para cnt en el mismo cluster
	        if (cluster[j] == i){ //Si el CNT pertenece al cluster
	        	for (int s=0;s<(n+1);s++){
	            	xmin = std::min((x[j*(n+1)+s]),xmin);
	            	xmax = std::max((x[j*(n+1)+s]),xmax);
	            }//end for s
	        }//end if
	    } //end for j
		if (((xmin<=0)&&(xmax>=Lx_strain))){ // si el cluster "atraviesa" el volumen representativo, el cluster ha percolado
			numero_clusters_percolados[mc] = numero_clusters_percolados[mc] + 1;
			indice_cluster_percolado_aux[(int)(numero_clusters_percolados[mc]-1)] = i;
		} //end if
	}//end for i			

	// DETERMINAR PARAMETROS DE INTERES EN CASO DE QUE HAYA PERCOLACION
	if((int)(numero_clusters_percolados[mc])>0){ //si algun cluster percola
		
		percolado[mc]=1;
		
		for (int i=0;i<(int)(numero_clusters_percolados[mc]);i++){ //se crea una lista con los cnt que pertenecen a algun cluster percolado
			for(int j=0;j<N;j++){
				if(cluster[j] == indice_cluster_percolado_aux[i]){ //si cnt j pertenece al cluster percolado i
					cnt_percolado[j] = 1;
				} //end if
			} //end for j
		} //end for i			
		
		//SE ELIMINA LA DISTANCIA ENTRE LOS CNT QUE NO ESTAN EN ALGUN CAMINO PERCOLADO 
		for (int k=0;k<numero_contactos;k++){   	
			int i = dmin[k*3];
			int j = dmin[k*3+1];
				if ( (cnt_percolado[i] == 0) || (cnt_percolado[j] == 0) )
					{dmin[k*3+2] = -1;}	
		}	
		
		//SE BUSCAN LOS NODOS EXTREMOS
		for (int i=0;i<N;i++){
			if (cnt_percolado[i] == 1){
				for (int s=0;s<(n+1);s++){
					if (x[i*(n+1)+s] < 0){
						nodo_extremo[i] = -1;
					}
					if (x[i*(n+1)+s] > Lx_strain){
						nodo_extremo[i] = 1;
					}
	            }//end for s
			} //end if percolado			
		} //end if for i 
//		cout<<"\nodo extremo:"<<endl;	imprimir1d(nodo_extremo,N);
		
		//SE ELIMINAN LAS DISTANCIAS DE LOS CNT QUE NO APORTAN A LA CONDUCTIVIDAD
		//(estos son los que no son extremos y estan en contacto con solo un cnt)
		eliminados = 1;
		while (eliminados>0){
			eliminados = 0;
			for (int i=0;i<N;i++){
				numero_contactos_cnt[i] = 0;
			}
			
			for (int k=0;k<numero_contactos;k++){			
				int i = dmin[k*3];
				int j = dmin[k*3+1];
				if ( (dmin[k*3+2]>0) && (dmin[k*3+2]<dt) ){
					numero_contactos_cnt[i] += 1;
					numero_contactos_cnt[j] += 1;	
				}
			}

			for (int k=0;k<numero_contactos;k++){			
				int i = dmin[k*3];
				int j = dmin[k*3+1];
				if ((dmin[k*3+2]>0) && (dmin[k*3+2]<dt)){
					if ( (nodo_extremo[i] == 0) && (numero_contactos_cnt[i] == 1) ){
						dmin[k*3+2] = d_max_inicial_strain;
						eliminados +=1;
					}
					if ( (nodo_extremo[j] == 0) && (numero_contactos_cnt[j] == 1) ){
						dmin[k*3+2] = d_max_inicial_strain;
						eliminados +=1;
					}
					
					if ( ((nodo_extremo[i] < 0) && (nodo_extremo[j] < 0)) 
					 ||  ((nodo_extremo[i] > 0) && (nodo_extremo[j] > 0)) ){ //caso de rama que parte desde nodo extremo 
						dmin[k*3+2] = d_max_inicial_strain;
						eliminados +=1;
					}
				}
			}
		} //end while
		//cout<<"\n dmin:"<<endl;	imprimir2d(dmin,numero_contactos,3);
		
		//SE LLENA LA MATRIX dmin_bb, LA CUAL CONTIENE LA DISTANCIA SOLO PARA LOS CNT QUE APORTAN A LA CONDUCTIVIDAD
		numero_edges_bb = 0;
		for (int k=0;k<numero_contactos;k++){	
			if ( (dmin[k*3+2]>0) && (dmin[k*3+2]<dt) ){
				dmin_bb[numero_edges_bb*3] = dmin[k*3];
				dmin_bb[numero_edges_bb*3+1] = dmin[k*3+1];
				dmin_bb[numero_edges_bb*3+2] = dmin[k*3+2];
				numero_edges_bb += 1;		
			}
		}
		// cout<<"\n dmin_bb:"<<endl;	imprimir2d(dmin_bb,numero_edges_bb,3);
		
		//SE DETERMINA LOS CNT QUE APORTAN A LA CONDUCTIVIDAD
		for (int k=0;k<numero_edges_bb;k++){			
			int i = dmin_bb[k*3];
			int j = dmin_bb[k*3+1];
			nodo_bb[i] = 1;
			nodo_bb[j] = 1;
		}
		
		//SE DETERMINA EL NUMERO DE CNT QUE APORTAN A LA CONDUCTIVIDAD, QUE EQUIVALE A LOS NODOS DE LA MATRIZ dmin_bb
		for (int i=0;i<N;i++){
			if (nodo_bb[i] > 0){
				numero_cnt_percolando[mc]++;
			}
		}
		
		numero_contactos_cnt_percolando[mc] = numero_edges_bb;	
		promedio_numero_contactos_cnt_percolando[mc] = 
			(double)((double)numero_contactos_cnt_percolando[mc]/(double)numero_cnt_percolando[mc]);	
		
		if (resistencias == 1)
		{	
			for (int i=0;i<N;i++)
			{
				if (numero_contactos_cnt[i] == 0)
				{
					nodo_extremo[i] = 0;
				} //end if
			} //end for 
			
			primer_nodo_max = -1;
			ipnm = 0;
			while (primer_nodo_max<=0){
				primer_nodo_max = nodo_extremo[ipnm];	
				ipnm++;
			}
			primer_nodo_max = ipnm-1;
			
			//Llenar la parte correspondiente a G
			indice_A_pre = 0;		
			for (int k=0;k<numero_edges_bb;k++){		
				int i = dmin_bb[k*3];
				int j = dmin_bb[k*3+1];
				double valor = G_tunel(dmin_bb[k*3+2],lambda);
				if (nodo_extremo[i]>0){ //unir nodos maximos
					i = primer_nodo_max;
				}
				if (nodo_extremo[j]>0){ //unir nodos maximos
					j = primer_nodo_max;
				}
				
				if (nodo_extremo[i]<0){ //considerar nodo minimo como tierra
					A_pre[indice_A_pre*3+0] = j;
					A_pre[indice_A_pre*3+1] = j;
					A_pre[indice_A_pre*3+2] = valor;
					indice_A_pre++;}
				else if (nodo_extremo[j]<0){ //considerar nodo minimo como tierra
					A_pre[indice_A_pre*3+0] = i;
					A_pre[indice_A_pre*3+1] = i;
					A_pre[indice_A_pre*3+2] = valor;
					indice_A_pre++;}
				else{
					A_pre[indice_A_pre*3+0] = i;
					A_pre[indice_A_pre*3+1] = i;
					A_pre[indice_A_pre*3+2] = valor;
					indice_A_pre++;
					A_pre[indice_A_pre*3+0] = j;
					A_pre[indice_A_pre*3+1] = j;
					A_pre[indice_A_pre*3+2] = valor;
					indice_A_pre++;
					A_pre[indice_A_pre*3+0] = i;
					A_pre[indice_A_pre*3+1] = j;
					A_pre[indice_A_pre*3+2] = -valor;
					indice_A_pre++;
					A_pre[indice_A_pre*3+0] = j;
					A_pre[indice_A_pre*3+1] = i;
					A_pre[indice_A_pre*3+2] = -valor;
					indice_A_pre++;
				} //end if
			}//end for
			
			//Llenar la parte correspondiente a B,C y D
			A_pre[indice_A_pre*3+0] = primer_nodo_max;
			A_pre[indice_A_pre*3+1] = (N+1)-1;
			A_pre[indice_A_pre*3+2] = 1;
			indice_A_pre++;
			
			A_pre[indice_A_pre*3+0] = (N+1)-1;
			A_pre[indice_A_pre*3+1] = primer_nodo_max;
			A_pre[indice_A_pre*3+2] = 1;
			indice_A_pre++;
			// cout<<"\n A_pre:"<<endl;	imprimir2dmin(A_pre,indice_A_pre,3);
				
			for (int k=0;k<indice_A_pre;k++){
				int i = A_pre[k*3+0];
				int j = A_pre[k*3+1];
				nodo_de_A[i] = 1;	
				nodo_de_A[j] = 1;	
			}
			
			numero_nodos_A = 0;
			for (int i=0;i<N+1;i++){
				if (nodo_de_A[i] == 1){
					indice_cnt_A[2*numero_nodos_A]=numero_nodos_A;
					indice_cnt_A[2*numero_nodos_A+1]=i;
					numero_nodos_A++;
				}
			}
			
			numero_elementos_A = 0;
			double* A=new double[numero_nodos_A*numero_nodos_A]();	
			double* b = new double[numero_nodos_A]();
			b[numero_nodos_A-1] = V;
			
			for (int k=0;k<indice_A_pre;k++){
				int i_pre = A_pre[k*3+0];
				int j_pre = A_pre[k*3+1];
				int encontrado_i = 0;
				int encontrado_j = 0;
				int i = 0;
				int j = 0;
				while (encontrado_i + encontrado_j < 2){
					if (indice_cnt_A[2*i+1] == i_pre){
						encontrado_i = 1;
					}
					else{
						i++;
					}
					
					if (indice_cnt_A[2*j+1] == j_pre){
						encontrado_j = 1;
					}
					else{
						j++;
					}
				}
				if (A[i*numero_nodos_A+j] == 0){
					numero_elementos_A++;
					}
				A[i*numero_nodos_A+j] = A[i*numero_nodos_A+j] + A_pre[k*3+2];
			}
			// cout<<"\n A: "<<endl;	imprimir2d(A,numero_nodos_A,numero_nodos_A);		
			
			int* Ai = new int[numero_elementos_A]();
			int* Ap = new int[(numero_nodos_A+1)]();
			double* Ax = new double[numero_elementos_A]();
			
			contador_Ap = 0;
			for (int col = 0;col<numero_nodos_A;col++){
				for (int row = 0;row<numero_nodos_A;row++){
					if(A[row*numero_nodos_A+col] != 0){
						Ai[contador_Ap] = row;
						Ax[contador_Ap] = A[row*numero_nodos_A+col];
						contador_Ap++;
					}
				Ap[col+1] = contador_Ap;
				}
			}
			
			for (int j=0; j<numero_nodos_A; j++)
			{
				bubble_sort(Ai,Ax,Ap[j],Ap[j+1]);
			}	

			klu_symbolic *Symbolic ;
			klu_numeric *Numeric ;
			klu_common Common ;
			
			klu_defaults (&Common) ;
			Symbolic = klu_analyze (numero_nodos_A, Ap, Ai, &Common) ;
			Numeric = klu_factor (Ap, Ai, Ax, Symbolic, &Common) ;
			klu_solve (Symbolic, Numeric, numero_nodos_A, 1, b, &Common) ;
			klu_free_symbolic (&Symbolic, &Common) ;
			klu_free_numeric (&Numeric, &Common) ;
			
			I = -(b[numero_nodos_A-1]);
			R[mc] = V/I;
			// cout << "R = "<<setprecision(5) << R[mc]<<endl; 

			delete[] Ap;		
			delete[] Ai;		
			delete[] Ax;		
			delete[] b;			
			delete[] A;
		} //end if resistencias
	} //end if clusters percolados
	
	if (guardar == 1){ //se guardan los datos en archivos .txt
		std::lock_guard<std::mutex> guard(myMutex);
		if (mc == 0){ //en la primera iteracion de mc se borran los datos anteriormente guardados
			std::ofstream x_file;
			x_file.open("resultados_resistencia_v4_0/xyz.txt",std::ios::out | std::ios::trunc | std::ios::binary);     
guardar_constantes(x_file,seed,N,L,n,n_aglomerados,n_aglomerados_x,n_aglomerados_y,n_aglomerados_z,p_aglomerado,theta_max,d0,dt,Lx,Ly,Lz,strain, poisson);
			x_file << "x=[";
			guardar2d(x_file,(double*)x,N,(n+1)); x_file << "];" <<"\n"<<endl;
			x_file << "y=[";
			guardar2d(x_file,(double*)y,N,(n+1)); x_file << "];" <<"\n"<< endl;
			x_file << "z=[";
			guardar2d(x_file,(double*)z,N,(n+1)); x_file << "];" <<"\n" << endl;
			x_file.close();
	
			std::ofstream xcm_file;
			xcm_file.open("resultados_resistencia_v4_0/xyz_cm.txt",std::ios::out | std::ios::trunc | std::ios::binary);     
guardar_constantes(xcm_file,seed,N,L,n,n_aglomerados,n_aglomerados_x,n_aglomerados_y,n_aglomerados_z,
p_aglomerado,theta_max, d0,dt,Lx,Ly,Lz,strain, poisson);
			xcm_file << "xcm=[";
			guardar2d(xcm_file,(double*)xcm,N,n); xcm_file << "];" <<"\n"<<endl;
			xcm_file << "ycm=[";
			guardar2d(xcm_file,(double*)ycm,N,n); xcm_file << "];" <<"\n"<<endl;
			xcm_file << "zcm=[";
			guardar2d(xcm_file,(double*)zcm,N,n); xcm_file << "];" <<"\n"<<endl;
			xcm_file.close();
			
			std::ofstream angulos_file;
angulos_file.open("resultados_resistencia_v4_0/angulos.txt",std::ios::out | std::ios::trunc | std::ios::binary);  
									
guardar_constantes(angulos_file,seed,N,L,n,n_aglomerados,n_aglomerados_x,n_aglomerados_y,n_aglomerados_z,
p_aglomerado,theta_max, d0,dt,Lx,Ly,Lz,strain, poisson);
			angulos_file << "theta=[";
			guardar2d(angulos_file,(double*)theta_suma,N,n); angulos_file << "];" <<"\n" <<endl;
			angulos_file << "phi=[";
			guardar2d(angulos_file,(double*)phi_suma,N,n); angulos_file << "];" <<"\n" << endl;
			angulos_file.close();
			
			std::ofstream cluster_file;
			cluster_file.open("resultados_resistencia_v4_0/resultados_cluster.txt",
std::ios::out | std::ios::trunc | std::ios::binary);    guardar_constantes(cluster_file,seed,N,L,n,n_aglomerados,n_aglomerados_x,n_aglomerados_y,n_aglomerados_z,
p_aglomerado,theta_max, d0,dt,Lx,Ly,Lz,strain, poisson);			
			guardar1d(cluster_file,cluster,N);
			guardar1d(cluster_file,indice_cluster_percolado_aux,numero_clusters_percolados[mc]);
			cluster_file.close();    
			
			std::ofstream graficar3d_matlab_file;
			std::stringstream stream;
			stream<<std::fixed<<"resultados_resistencia_v4_0/s"
<<std::setprecision(0)<<seed<<" "<<N <<" "<<L<<" "<<n<<" "<<n_aglomerados<<" "<<n_aglomerados_x
<<" "<<n_aglomerados_y<<" "<<n_aglomerados_z<<" "
			<<std::setprecision(2)<<p_aglomerado<<" "
			<<std::setprecision(2) <<theta_max<<" "
			<<std::setprecision(0)<<d0<<" "<<dvdw<<" "<<dt<<" "<<Lx/L<<" "<<Ly/L<<" "<<Lz/L<<" "
			<<std::setprecision(2)<<strain<<" "<<poisson<<"&m";
			std::string file_name = stream.str();
			std::replace( file_name.begin(), file_name.end(), '.', '_');
			std::replace( file_name.begin(), file_name.end(), ' ', '_');
			std::replace( file_name.begin(), file_name.end(), '&', '.');
			
			graficar3d_matlab_file.open(file_name,std::ios::out | std::ios::trunc | std::ios::binary); 
graficar3d_matlab_file << "constantes=[";  guardar_constantes(graficar3d_matlab_file,seed,N,L,n,n_aglomerados,n_aglomerados_x,n_aglomerados_y,n_aglomerados_z,p_aglomerado, theta_max,d0,dt,Lx,Ly,Lz,strain, poisson);
			graficar3d_matlab_file << "];" <<"\n"<<endl;					
			graficar3d_matlab_file << "x=[";
			guardar2d(graficar3d_matlab_file,(double*)x,N,(n+1)); graficar3d_matlab_file << "];" <<"\n"<<endl;
			graficar3d_matlab_file << "y=[";
			guardar2d(graficar3d_matlab_file,(double*)y,N,(n+1)); graficar3d_matlab_file << "];" <<"\n"<< endl;
			graficar3d_matlab_file << "z=[";
			guardar2d(graficar3d_matlab_file,(double*)z,N,(n+1)); graficar3d_matlab_file << "];" <<"\n" << endl;
			graficar3d_matlab_file << "grafico3d(x,y,z,constantes);" << endl;
			graficar3d_matlab_file << "histogramas(x,y,z,constantes);" << endl;
			graficar3d_matlab_file.close();  			

			std::ofstream spice_file;
			spice_file.open("resultados_resistencia_v4_0/spice.cir",std::ios::out | std::ios::trunc | std::ios::binary);     
			spice(spice_file,(double*)dmin_bb,numero_edges_bb,(int*)nodo_extremo,N,lambda);
			spice_file.close();			
		} //end if mc=0
		else{ //para las siguientes iteraciones de mc se agregan los datos al final
			std::ofstream x_file;
			x_file.open("resultados_resistencia_v4_0/xyz.txt",std::ios::out | std::ios::app | std::ios::binary);    
guardar_constantes(x_file,seed,N,L,n,n_aglomerados,n_aglomerados_x,n_aglomerados_y,n_aglomerados_z,
p_aglomerado, theta_max,d0,dt,Lx,Ly,Lz,strain, poisson);
			x_file << "x=[";
			guardar2d(x_file,(double*)x,N,(n+1)); x_file << "];" <<"\n"<<endl;
			x_file << "y=[";
			guardar2d(x_file,(double*)y,N,(n+1)); x_file << "];" <<"\n"<< endl;
			x_file << "z=[";
			guardar2d(x_file,(double*)z,N,(n+1)); x_file << "];" <<"\n" << endl;
			x_file.close();
	
			std::ofstream xcm_file;
			xcm_file.open("resultados_resistencia_v4_0/xyz_cm.txt",std::ios::out | std::ios::app | std::ios::binary);   
guardar_constantes(xcm_file,seed,N,L,n,n_aglomerados,n_aglomerados_x,n_aglomerados_y,n_aglomerados_z,
p_aglomerado, theta_max,d0,dt,Lx,Ly,Lz,strain, poisson);
			xcm_file << "xcm=[";
			guardar2d(xcm_file,(double*)xcm,N,n); xcm_file << "];" <<"\n"<<endl;
			xcm_file << "ycm=[";
			guardar2d(xcm_file,(double*)ycm,N,n); xcm_file << "];" <<"\n"<<endl;
			xcm_file << "zcm=[";
			guardar2d(xcm_file,(double*)zcm,N,n); xcm_file << "];" <<"\n"<<endl;
			xcm_file.close();
			
			std::ofstream angulos_file;
			angulos_file.open("resultados_resistencia_v4_0/angulos.txt",
std::ios::out | std::ios::app | std::ios::binary);  
guardar_constantes(angulos_file,seed,N,L,n,n_aglomerados,n_aglomerados_x,n_aglomerados_y,n_aglomerados_z,
p_aglomerado, theta_max,d0,dt,Lx,Ly,Lz,strain, poisson);
			angulos_file << "theta=[";
			guardar2d(angulos_file,(double*)theta_suma,N,n); angulos_file << "];" <<"\n" <<endl;
			angulos_file << "phi=[";
			guardar2d(angulos_file,(double*)phi_suma,N,n); angulos_file << "];" <<"\n" << endl;
			angulos_file.close();
			
			std::ofstream cluster_file;
			cluster_file.open("resultados_resistencia_v4_0/resultados_cluster.txt",
std::ios::out | std::ios::app | std::ios::binary);  
guardar_constantes(cluster_file,seed,N,L,n,n_aglomerados,n_aglomerados_x,n_aglomerados_y,n_aglomerados_z,
p_aglomerado, theta_max,d0,dt,Lx,Ly,Lz,strain, poisson);
			guardar1d(cluster_file,cluster,N);
			guardar1d(cluster_file,indice_cluster_percolado_aux,numero_clusters_percolados[mc]);
			cluster_file.close();  	
		} // end else
	} // end if guardar
	
	double promedio_conexiones_A;
	if (numero_nodos_A == 0)
		{promedio_conexiones_A = 0;} 
	else
		{promedio_conexiones_A = (double)((double)numero_elementos_A/(double)numero_nodos_A);}
			
	//Se guardan los resultados para cada seed
	if (guardar == 2){
		std::ofstream resultados_por_seed_file;
		resultados_por_seed_file.open("resultados_resistencia_v4_0/resultados_por_seed.txt",
std::ios::out | std::ios::app | std::ios::binary);     
		resultados_por_seed_file<< std::fixed
			<<setprecision(0)<< seed <<" "<< N <<" "<< L <<" "<< n <<" "
			<<setprecision(0)<< n_aglomerados <<" "<< n_aglomerados_x <<" "<< n_aglomerados_y <<" "<< n_aglomerados_z <<" "
			<<setprecision(2)<< p_aglomerado <<" "<< theta_max <<" "
			<<setprecision(0)<< d0 <<" "<< dvdw <<" "<< dt <<" "<< Lx <<" "<< Ly <<" "<< Lz <<" "
			<<setprecision(2)<< strain <<" "<< poisson <<" "
			<<setprecision(0)<<R_pol<<" "
			<<setprecision(0)<<percolado[mc]<<" "
			<<setprecision(0)<<numero_cnt_percolando[mc]<<" "
			<<setprecision(0)<<numero_contactos_cnt_percolando[mc]<<" "
			<<setprecision(1)<<promedio_numero_contactos_cnt_percolando[mc]<<" "
			<<setprecision(1)<<((double)(numero_contactos*100/N)/N) <<" "
			<<setprecision(0)<<numero_nodos_A<<" "
			<<setprecision(0)<<numero_elementos_A<<" "
			<<setprecision(1)<<promedio_conexiones_A<<" "
			<<setprecision(6)<< R[mc]<<" " 
			<< std::endl;    
		resultados_por_seed_file.close(); 
	}
	
	if (imprimir==2)
	{
		if (mc % 50 == 0)
		{std::cout <<std::fixed
			<<setprecision(0)<<"sed:" << seed 
			<<setprecision(0)<<" Box:"<< Lx/L 
			<<setprecision(0)<<" n:"<<n
			<<setprecision(0)<<" N:" << N 
			<<setprecision(0)<<" n_a:" << n_aglomerados 
			<<setprecision(2)<<" p_a:" << p_aglomerado
			<<setprecision(2)<<" str:" << strain 
			<<setprecision(0)<<" cnt_p:"<<numero_cnt_percolando[mc]<<" "
			<<setprecision(0)<<" n_cont_p:"<<numero_contactos_cnt_percolando[mc]<<" "
			<<setprecision(1)<<" prm_cont_p:"<<promedio_numero_contactos_cnt_percolando[mc]<<" "
			<<setprecision(2)<<" %cont:"<<((double)(numero_contactos*100/N)/N) <<" "
			<<setprecision(0)<<" n_nodo_A:"<<numero_nodos_A<<" "
			<<setprecision(0)<<" n_elm_A:"<<numero_elementos_A<<" "
			<<setprecision(1)<<" prm_cnx_A:"<<promedio_conexiones_A<<" "
			<<setprecision(6)<<" R:"<< R[mc] 
			<<std::endl; 
		}
	} //end if imprimir
	
	//liberar memoria de los array creados	
	delete[] x;	
	delete[] y;	
	delete[] z;
	delete[] theta_suma;	
	delete[] phi_suma;	
	delete[] xcm;	
	delete[] ycm;	
	delete[] zcm;
	delete[] kx;	
	delete[] ky;	
	delete[] kz;
	delete[] rcm_i;	
	delete[] rcm_j;		
	delete[] dcmkm;
	delete[] rki;
	delete[] rmj;
	delete[] ri;
	delete[] rj;
	delete[] p_i;
	delete[] pj;
	delete[] rij;
	delete[] rji;
	delete[] v;	
	delete[] cluster;
	delete[] indice_cluster_percolado_aux;
	delete[] cnt_percolado; 
	delete[] dmin;	
	delete[] numero_contactos_cnt;
	delete[] nodo_extremo;
	delete[] dmin_bb;
	delete[] nodo_bb;	
	delete[] A_pre;
	delete[] nodo_de_A;	
	delete[] indice_cnt_A;	
} // end funcion iteracion_mc