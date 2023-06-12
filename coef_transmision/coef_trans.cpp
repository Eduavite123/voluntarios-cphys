#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <stdlib.h>

#define MAX 1000
#define PI 3.14159265

using namespace std;

int main(){
    //------------------------DECLARACION DE VARIABLES-------------------------------------
    complex<double> phi[MAX]; //Solución ec. Schrodinger  
    double norma; //Norma de la f. de onda
    int j,n; //posición y tiempo
    int t; //Tiempo máximo
    //Asumimos los parámetros de discretización espacial y temporal h=1 y s=1, respectivamente

    //Valores iniciales:
    int N; //tamaño pozo y nº ciclos (nc=1,...,N/4 para que haya 4 ciclos)
    float nc;
    float lambda; //proporcional a la altura del pozo (p.e. lambda=0.3)
    double x_0, sigma; //centro y anchura de la gaussiana de la f. onda inicial
    //Valores iniciales a generar:
    double S; //Espaciado temporal reescrito
    double k_0; //Proporcional a la altura del potencial 
    double V[MAX]; //Potencial según posición
    //Parámetros para la resolución del sistema de ecuaciones
    complex<double> alpha[MAX], A_neg, A_0[MAX], A_pos, gamma[MAX], beta[MAX], b[MAX], xi[MAX];
    complex<double> i (0.0,1.0); 

    //Cálculo del coeficiente de transmisión
    double mT=0; //nº de veces que se encuentra la partícula a la derecha del potencial
    double P_D=0; //Probabilidad de encontrar la partícula a la derecha de la barrera
    int n_D=0; //tiempo n en el que se encuentra el primer máximo local de P_D(t)
    double prev=0, PD_nD;
    double p, K;
    //-------------------------------------------------------------------------------------

    //--------------------------CONDICIONES INICIALES--------------------------------------
    t=MAX*2;
    //Parámetros del pozo:
    N=MAX; //Tamaño del pozo  
    nc=N/15.0; //# ciclos. Lo elegimos así para que haya 4 ciclos

    //Parámetros del potencial:
    lambda=1; //Altura del potencial
    k_0=2*PI*nc/(N*1.0);
    //Definimos el potencial
    for (j = 0; j < N; j++)
    {
        if (j<3.0/5.0*(double)N && j>2.0/5.0*(double)N)
        {
            V[j]=lambda*k_0*k_0;
        }
        else V[j]=0;        
    }

    //Condiciones de contorno:
    phi[0]=0.0;
    phi[N]=0.0;
    //Asumimos la función de onda inicial como una gaussiana
    /*Primero decidimos el centro y la anchura de la gaussiana 
    (conviene que la anchura sea mucho menor que el pozo y que esté centrada entre la primera 
    pared y el potencial)*/
    x_0=(double)N/4.0;
    sigma=(double)N/16.0;
    norma=0.0;
    for (j = 1; j < N-1; j++)
    {
        phi[j]=exp(i*k_0*(double)j)*exp((-((double)j-x_0)*((double)j-x_0))/(2*sigma*sigma));
        norma=norma+abs(phi[j])*abs(phi[j]);
    }
    //Normalizamos la f. de onda para que después salga la norma cercana a 1
    for(j=0;j<N;j++)
    {
        phi[j]=phi[j]/sqrt(norma);
    }

    //Parámetros para la resolución del sistema de ecuaciones:
    S=1.0/(4.0*k_0*k_0); //parámetro de espaciado temporal
    A_neg=A_pos=1;
    for(j=0;j<N;j++)
    {
        A_0[j]=-2.0+2.0*i/S-V[j];
    }
    alpha[N-1]=0.0;
    for(j=0;j<=N-2;j++)
    {
        alpha[(N-2)-j]=-1.0/(A_0[(N-1)-j]+alpha[(N-1)-j]);
    }
    //-------------------------------------------------------------------------------------

    //--------------------------ALGORITMO EC. SCHRODINGER----------------------------------
    for(n=0;n<t;n++)
    {
        P_D=0;
        for (j = N/5; j < N; j++)
        {
            P_D=P_D+abs(phi[j])*abs(phi[j]);
        }
        if (prev>P_D)
        {
            n_D=n-1;
            PD_nD=prev;
        }
        else
        {
            prev=P_D;
        }      
        
        //Cálculo de b
        for(j=0;j<N;j++)
        {
            b[j]=4.0*i*phi[j]/S;
        }
        //Cálculo de beta
        beta[N-1]=0; 
        for(j=0;j<=N-2;j++)
        {
            beta[(N-2)-j]=(b[(N-1)-j]-beta[(N-1)-j])/(A_0[(N-1)-j]+alpha[(N-1)-j]);
        }

        //Cálculo de Xi
        xi[0]=xi[N-1]=0;
        for(j=0;j<N-1;j++)
        {
            xi[j+1]=alpha[j]*xi[j]+beta[j];
        }

        //Cálculo de phi en el siguiente instante de tiempo y la norma
        for(j=0;j<N;j++)
        {
            phi[j]=xi[j]-phi[j];
        }
    }

    cout << PD_nD << "  " << n_D << endl;

    srand(time(NULL));
    for (int i = 0; i < 1000; i++)
    {
        p=(double)rand()/RAND_MAX;
        if (p>PD_nD)
        {
            mT=mT+1;
        }
        else if (p<PD_nD)
        {
            mT=mT;
        }
        
    }

    K=mT/1000.0;
    cout << mT << endl;
    cout << "K= " << K << endl;

    return 0;
}