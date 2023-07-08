/*Con este programa probaremos la eficacia del superordenador simulando simultáneamente varios sistemas
con condiciones iniciales ligeramente distintas y probar nuevamente la sensibilidad de este
tipo de sistemas*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono> //calcularemos el tiempo de ejecución

using namespace std;
using namespace std::chrono;
#define g 9.81
#define PI 3.14159265

double f0(double O1, double O2, double p_o1, double p_o2);
double f1(double O1, double O2, double p_o1, double p_o2);
double f2(double O1, double O2, double p_o1, double p_o2);
double f3(double O1, double O2, double p_o1, double p_o2);

int main(){
    //------------------------DECLARACION DE VARIABLES-------------------------------------
    double y[4]; //funciones a resolver de las ecuaciones diferenciales
    //   y_0 = \theta_1     y_1 = \theta_2     y_2 = p_{\theta_1}  y_3 = p_{\theta_2}
    float t, h, tmax, iter; //Tiempo y pasos temporales | Tiempo máx y nº iteraciones
    double k1[4], k2[4], k3[4], k4[4];
    double aux[4];
    double dO1, dO2; //Velocidades angulares
    double H; //Hamiltoniano
    //-------------------------------------------------------------------------------------

    //Paso de tiempo y tiempo máximo
    h=0.1;
    tmax=1000;
    iter=tmax/h;

    auto start = high_resolution_clock::now();

    //SISTEMA 1
    //--------------------------CONDICIONES INICIALES--------------------------------------
    //Posiciones iniciales de las masas
    y[0]=PI/2.0;
    y[1]=PI/2.0;
    //Momentos iniciales
    y[2]=0.0;
    y[3]=0.0;
    //-------------------------------------------------------------------------------------

    //-------------------------------RUNGE_KUTTA-------------------------------------------
    ofstream pos1;
    pos1.open("joel1.dat");
    for (int i = 0; i < iter; i++)
    {
        
        //Escribimos los datos en el archivo
        //Simulación péndulo doble
        pos1 << sin(y[0]) << ", " << -cos(y[0]) << endl;
        pos1 << sin(y[0])+sin(y[1]) << ", " << -cos(y[0])-cos(y[1]) << endl;
        pos1 << endl; 

        //Evaluamos k1
        k1[0]=h*f0(y[0], y[1], y[2], y[3]);
        k1[1]=h*f1(y[0], y[1], y[2], y[3]);
        k1[2]=h*f2(y[0], y[1], y[2], y[3]);
        k1[3]=h*f3(y[0], y[1], y[2], y[3]);

        //Un vector auxiliar facilitará el cálculo
        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k1[j]/2.0;
        }

        //Evaluamos k2
        k2[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k2[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k2[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k2[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);

        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k2[j]/2.0;
        }

        //Evaluamos k3
        k3[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k3[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k3[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k3[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);        
        
        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k3[j];
        }

        //Evaluamos k4
        k4[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k4[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k4[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k4[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);   

        //Finalmente evaluamos y_n(t+h)
        for (int k = 0; k < 4; k++)
        {
            y[k]=y[k]+1.0/6.0*(k1[k]+2*k2[k]+2*k3[k]+k4[k]);
        }

        t=t+h; 
  
    }
    pos1.close();
    //-------------------------------------------------------------------------------------

    //SISTEMA 2
        //--------------------------CONDICIONES INICIALES--------------------------------------
    //Posiciones iniciales de las masas
    y[0]=PI/2.0;
    y[1]=PI/2.001;
    //Momentos iniciales
    y[2]=0.0;
    y[3]=0.0;
    //-------------------------------------------------------------------------------------

    //-------------------------------RUNGE_KUTTA-------------------------------------------
    ofstream pos2;
    pos2.open("joel2.dat");
    for (int i = 0; i < iter; i++)
    {
        
        //Escribimos los datos en el archivo
        //Simulación péndulo doble
        pos2 << sin(y[0]) << ", " << -cos(y[0]) << endl;
        pos2 << sin(y[0])+sin(y[1]) << ", " << -cos(y[0])-cos(y[1]) << endl;
        pos2 << endl; 

        //Evaluamos k1
        k1[0]=h*f0(y[0], y[1], y[2], y[3]);
        k1[1]=h*f1(y[0], y[1], y[2], y[3]);
        k1[2]=h*f2(y[0], y[1], y[2], y[3]);
        k1[3]=h*f3(y[0], y[1], y[2], y[3]);

        //Un vector auxiliar facilitará el cálculo
        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k1[j]/2.0;
        }

        //Evaluamos k2
        k2[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k2[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k2[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k2[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);

        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k2[j]/2.0;
        }

        //Evaluamos k3
        k3[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k3[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k3[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k3[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);        
        
        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k3[j];
        }

        //Evaluamos k4
        k4[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k4[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k4[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k4[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);   

        //Finalmente evaluamos y_n(t+h)
        for (int k = 0; k < 4; k++)
        {
            y[k]=y[k]+1.0/6.0*(k1[k]+2*k2[k]+2*k3[k]+k4[k]);
        }

        t=t+h; 
  
    }
    pos2.close();
    //-------------------------------------------------------------------------------------

    //SISTEMA 3
        //--------------------------CONDICIONES INICIALES--------------------------------------
    //Posiciones iniciales de las masas
    y[0]=PI/2.0;
    y[1]=PI/2.0005;
    //Momentos iniciales
    y[2]=0.0;
    y[3]=0.0;
    //-------------------------------------------------------------------------------------

    //-------------------------------RUNGE_KUTTA-------------------------------------------
    ofstream pos3;
    pos3.open("joel3.dat");
    for (int i = 0; i < iter; i++)
    {
        
        //Escribimos los datos en el archivo
        //Simulación péndulo doble
        pos3 << sin(y[0]) << ", " << -cos(y[0]) << endl;
        pos3 << sin(y[0])+sin(y[1]) << ", " << -cos(y[0])-cos(y[1]) << endl;
        pos3 << endl; 

        //Evaluamos k1
        k1[0]=h*f0(y[0], y[1], y[2], y[3]);
        k1[1]=h*f1(y[0], y[1], y[2], y[3]);
        k1[2]=h*f2(y[0], y[1], y[2], y[3]);
        k1[3]=h*f3(y[0], y[1], y[2], y[3]);

        //Un vector auxiliar facilitará el cálculo
        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k1[j]/2.0;
        }

        //Evaluamos k2
        k2[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k2[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k2[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k2[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);

        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k2[j]/2.0;
        }

        //Evaluamos k3
        k3[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k3[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k3[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k3[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);        
        
        for (int j = 0; j < 4; j++)
        {
            aux[j]=y[j]+k3[j];
        }

        //Evaluamos k4
        k4[0]=h*f0(aux[0],aux[1],aux[2],aux[3]);
        k4[1]=h*f1(aux[0],aux[1],aux[2],aux[3]);
        k4[2]=h*f2(aux[0],aux[1],aux[2],aux[3]);
        k4[3]=h*f3(aux[0],aux[1],aux[2],aux[3]);   

        //Finalmente evaluamos y_n(t+h)
        for (int k = 0; k < 4; k++)
        {
            y[k]=y[k]+1.0/6.0*(k1[k]+2*k2[k]+2*k3[k]+k4[k]);
        }

        t=t+h; 
  
    }
    pos3.close();
    //-------------------------------------------------------------------------------------

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
 
    cout << "Tiempo de ejecución: "<< duration.count() << " microsegundos" << endl;

    return 0;
}

double f0(double O1, double O2, double p_o1, double p_o2){
    double f0;
    f0=(p_o1-p_o2*cos(O1-O2))/(1+pow(sin(O1-O2),2));
    return f0;
}

double f1(double O1, double O2, double p_o1, double p_o2){
    double f1;
    f1=(-p_o1*cos(O1-O2)+2*p_o2)/(1+pow(sin(O1-O2),2));
    return f1;
}

double f2(double O1, double O2, double p_o1, double p_o2){
    double f2, h1, h2;
    h1=(p_o1*p_o2*sin(O1-O2))/(1+pow(sin(O1-O2),2));
    h2=(pow(p_o1,2)+2*pow(p_o2,2)-2*p_o1*p_o2*cos(O1-O2))/(2*pow(1+pow(sin(O1-O2),2),2));
    f2=-2*g*sin(O1)-h1+h2*sin(2*(O1-O2));
    return f2;
}

double f3(double O1, double O2, double p_o1, double p_o2){
    double f3, h1, h2;
    h1=(p_o1*p_o2*sin(O1-O2))/(1+pow(sin(O1-O2),2));
    h2=(pow(p_o1,2)+2*pow(p_o2,2)-2*p_o1*p_o2*cos(O1-O2))/(2*pow(1+pow(sin(O1-O2),2),2));
    f3=-g*sin(O2)+h1-h2*sin(2*(O1-O2));
    return f3;
}