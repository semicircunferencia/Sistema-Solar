/*
ALGORITMO DE VERLET
Desarrollamos r_i(t+h) y v_i(t+h) en serie de Taylor. Todo son magnitudes vectoriales. Queda:

r_i(t+h) = r_i(t) + h*v_i(t) + h^2*a_i(t)/2
v_i(t+h) = v_i(t) + h*(a_i(t) + a_i(t+h))/2

Los pasos a seguir son:
1) (Solo la primera vez, ver paso 3) Calcular a_i(t) a partir de los r_i(t), usando la ley de gravitación universal.
2) Calcular r_i(t+h) a partir de r_i(t), v_i(t), a_i(t).
3) Calcular a_i(t+h) con la ley de gravitación. ¡ Tener en cuenta que ya no hace falta calcularlo en la siguiente iteración !
4) Calcular v_i(t+h) usando los datos anteriores.


REESCALAMOS LAS VARIABLES
Escribimos r en UA, m en masas solares, y reescalamos el tiempo como t'=(GM_s/d^3)^(1/2)t, donde M_s es la masa del sol y
d= 1 UA. En estas unidades, G=1.

CONDICIONES INICIALES
Recuerda que podemos restringir el sistema a dos dimensiones. Cada r_i, v_i tendrá dos componentes. Como
condiciones iniciales, podemos colocar todos los planetas en el eje x y todas las velocidades en el y (perpendiculares).

OBJETIVOS:
1) Comprobar que se obtienen órbitas elípticas.
2) Comprobar que el periodo T de cada planeta se asemeja al real.
3) Observar si es un sistema estable, perturbando ligeramente cada planeta.
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

// Número de cuerpos
#define N 9

// Step temporal
#define h 1e-3

// Número de iteraciones
#define iter 1.1e6

// Constantes para renormalizar los parámetros
#define ctetiempo 0.1593 // Si multiplicas t por esto, sale el tiempo en años
using namespace std;

// Declaro las funciones
void leercondiniciales(string nombre, double masas[], double posiciones[][2], double velocidades[][2]);
void calcularacelentmash(double posiciones[][2], double velocidades[][2], double acelent[][2], double acelentmash[][2],
double masas[]);
void iteracionVerlet(double posiciones[][2], double velocidades[][2], double acelent[][2], double acelentmash[][2],
double masas[]);
void Guniv(double posiciones[][2], double aceleraciones[][2], double masas[]);
double energiakin(double velocidades[][2], double masas[]);
double energiapot(double posiciones[][2], double masas[]);
double momentoangular(double posiciones[][2], double velocidades[][2], double masas[]);
void integralangulo(double angulos[], double posiciones[][2], double velocidades[][2], bool vueltadada[]);
double rVerlet(double r, double v, double a);
double vVerlet(double v, double a, double a2);


/**************************************************** FUNCIÓN PRINCIPAL ****************************************************/
int main(void) {

    // Vectores de masas, posiciones, velocidades, aceleraciones en t y en t+h
    double masas[N];
    double posiciones[N][2];
    double velocidades[N][2];
    double acelent[N][2];
    double acelentmash[N][2];

    // Para controlar los periodos
    double angulos[N];
    bool vueltadada[N];
    bool periodoimpreso[N];

    // Leo las condiciones iniciales
    leercondiniciales("Condiniciales.txt", masas, posiciones, velocidades);

    // Inicializo los ángulos a 0 y los booleanos a falso
    for(int i=1; i<N; i++) angulos[i]=0;
    for(int i=1; i<N; i++) vueltadada[i]=false;
    for(int i=1; i<N; i++) periodoimpreso[i]=false;
    
    // Calculo las aceleraciones en el instante 0 y 1. Única vez que llamo a estas funciones aqui
    Guniv(posiciones, acelent, masas);
    calcularacelentmash(posiciones, velocidades, acelent, acelentmash, masas);

    // Abro los ficheros: uno para guardar, uno para el vídeo, uno para la energía y momento y otro para los periodos
    ofstream datos;
    ofstream datospython;
    ofstream datosenergiaymomento;
    ofstream datosperiodos;

    datos.open("Todo.dat");
    datospython.open("planets_data.dat");
    datosenergiaymomento.open("Energia_y_momento.dat");
    datosperiodos.open("Periodos.dat");

    // Número de iteraciones en el tiempo
    for(int j=0; j<=iter; j++) {
        // Para cada planeta, pego los datos en los ficheros, pero solo para cada 100 iteraciones
        if(j%100==0) {
            // Escribo el tiempo en el fichero con todo
            datos << ctetiempo*j*h << "\n";
            for(int i=0; i<N; i++) {
                // El fichero con todo
                for(int k=0; k<2; k++) datos << setw(15) << posiciones[i][k];
                for(int k=0; k<2; k++) datos << setw(15) << velocidades[i][k];
                for(int k=0; k<2; k++) datos << setw(15) << acelent[i][k];
        
                datos << "\n";

                // El fichero de Python
                datospython << posiciones[i][0]-posiciones[3][0] << "," << posiciones[i][1]-posiciones[3][1] << "\n";

            }

            datos << "\n";
            datospython << "\n";
            
            // Escribo los datos de energía y momento angular en su fichero
            datosenergiaymomento << setw(15) << j*h << setw(15) <<
            energiakin(velocidades, masas)+energiapot(posiciones, masas) << setw(15) <<
            momentoangular(posiciones, velocidades, masas) << "\n";

        }

        // Añado la contribución al ángulo y compruebo si se ha completado una vuelta
        integralangulo(angulos, posiciones, velocidades, vueltadada);

        // Imprimo si la vuelta está dada y el periodo no está impreso todavía
        for(int i=1; i<N; i++) {
            if(vueltadada[i] && !periodoimpreso[i]) {
                datosperiodos << "Planeta " << i << setw(15) << ctetiempo*j*h << "\n";

                // El periodo está impreso
                periodoimpreso[i]=true;
            }
        }
        
        

        // Calculo los nuevos parámetros
        iteracionVerlet(posiciones, velocidades, acelent, acelentmash, masas);
    }

    datos.close();
    datospython.close();
    datosenergiaymomento.close();
    datosperiodos.close();

    return 0;
}
/***************************************************************************************************************************/


/*Función leercondiniciales. Lee el fichero de condiciones iniciales, de tres columnas, una con las masas, otra 
con las posiciones y otra con las velocidades. Lo asigna a los vectores o arrays correspondientes*/
void leercondiniciales(string nombre, double masas[], double posiciones[][2], double velocidades[][2]) {

    // Defino radio y ángulo
    double radio[N];
    double angulo[N]; 

    // Intento abrir el fichero
    ifstream fichero;
    fichero.open(nombre);

    // Si encuentra el archivo, copia su contenido en los tres vectores
    if(fichero.is_open()) {
        for(int j=0; j<N; j++) {
            // Lee la masa
            fichero >> masas[j];

            // Lee la posición inicial a partir del radio y ángulo
            fichero >> radio[j];
            fichero >> angulo[j];

            posiciones[j][0]=radio[j]*cos(angulo[j]);
            posiciones[j][1]=radio[j]*sin(angulo[j]);

            // Lee las velocidades iniciales en los ejes x e y, respectivamente
            fichero >> velocidades[j][0];
            fichero >> velocidades[j][1];

        }

        // Para mantener el sistema centrado, calculo la velocidad inicial del Sol de manera que la del CM del sistema sea 0
        // Calculo el momento total de los planetas
        double momento[2];
        for(int k=0; k<2; k++) momento[k]=0;
        for(int j=1; j<N; j++) {
            // Componente k
            for(int k=0; k<2; k++) momento[k]+=-masas[j]*velocidades[j][k];
        }

        for(int k=0; k<2; k++) velocidades[0][k]=momento[k]/masas[0];


        fichero.close();
    }

    // Si no, devuelve un mensaje de error
    else {
        cout << "Error. No se pudo abrir el archivo.\n";

    }

    return;
}


/*Función calcularacelentmash. Calcula las posiciones, velocidades y aceleraciones en el instante h
a partir de las anteriores, en 0*/
void calcularacelentmash(double posiciones[][2], double velocidades[][2], double acelent[][2], double acelentmash[][2],
double masas[]) {

    double posicionesenh[N][2];

    // Calculo la posición en el instante h para cada planeta
    for(int i=0; i<N; i++) {
        // Componentes x e y
        for(int k=0; k<2; k++) {
            posicionesenh[i][k]=rVerlet(posiciones[i][k], velocidades[i][k], acelent[i][k]);
        }

    }

    // A partir de la ley de gravitación universal, calculo la nueva aceleración en t+h
    Guniv(posicionesenh, acelentmash, masas);

    return;
}

/*Función iteracionVerlet. Calcula las posiciones, velocidades y aceleraciones en el instante t+h
a partir de las anteriores, en t*/
void iteracionVerlet(double posiciones[][2], double velocidades[][2], double acelent[][2], double acelentmash[][2],
double masas[]) {

    // Calculo la posición en el instante t+h para cada planeta
    for(int i=0; i<N; i++) {
        // Componentes x e y
        for(int k=0; k<2; k++) {
            posiciones[i][k]=rVerlet(posiciones[i][k], velocidades[i][k], acelent[i][k]);
        }

    }

    // Aceleración en t+h ahora es en t
    for(int i=0; i<N; i++) {
        // Componentes x e y
        for(int k=0; k<2; k++) {
            acelent[i][k]=acelentmash[i][k];
        }

    }

    // A partir de la ley de gravitación universal, calculo la nueva aceleración en t+h
    Guniv(posiciones, acelentmash, masas);

    // Calculo la velocidad en el instante t+h
    for(int i=0; i<N; i++) {
        // Componentes x e y
        for(int k=0; k<2; k++) {
            velocidades[i][k]=vVerlet(velocidades[i][k], acelent[i][k], acelentmash[i][k]);
        }

    }

    return;
}

/*Función Guniv. Calcula ambas componentes de la aceleración de todos los planetas a partir de las masas y las posiciones.*/
void Guniv(double posiciones[][2], double aceleraciones[][2], double masas[]) {

    for(int i=0; i<N; i++) {

        // Inicializo las componentes de las aceleraciones a[i][0] y a[i][1]
        aceleraciones[i][0]=0; aceleraciones[i][1]=0;


        // Sumo la contribución de la masa j, j distinto de i
        for (int j=0; j<N; j++) {
            if(j==i) continue;
            else {
                // k=1 es la componente x, k=2 es la componente y
                for(int k=0; k<2; k++) {
                    aceleraciones[i][k]+=-masas[j]*
                    (posiciones[i][k]-posiciones[j][k])/
                    pow(pow(posiciones[i][0]-posiciones[j][0],2)+pow(posiciones[i][1]-posiciones[j][1],2),1.5);
                }
            }
        
        }

    }
    return;
}

/*Función energiakin. Calcula la energía cinética del sistema*/
double energiakin(double velocidades[][2], double masas[]) {
    double energia;

    energia=0;
    for(int i=0; i<N; i++) {

        // Sumo la energía cinética de la masa i
        energia+=masas[i]*(velocidades[i][0]*velocidades[i][0]+velocidades[i][1]*velocidades[i][1])/2;
    }

    return energia;
    
}

/*Función energiapot. Calcula la energía potencial del sistema en un conjunto de posiciones*/
double energiapot(double posiciones[][2], double masas[]) {
    double energia;

    energia=0;
    for(int i=1; i<N; i++) {

        // Sumo la interacción de la masa k<i con la masa i
        for(int k=0; k<i; k++) energia+=-masas[i]*masas[k]/sqrt((posiciones[k][0]-posiciones[i][0])*
        (posiciones[k][0]-posiciones[i][0])+(posiciones[k][1]-posiciones[i][1])*(posiciones[k][1]-posiciones[i][1]));
    }

    return energia;
    
}

/*Función momentoangular. Calcula el momento angular a partir de las masas, posiciones y velocidades*/
double momentoangular(double posiciones[][2], double velocidades[][2], double masas[]) {
    double momento;

    momento=0;
    for(int i=0; i<N; i++) {

        // Sumo la contribución de la masa i. Sabemos que L=m*(x*v_y - y*v_x)
        momento+=masas[i]*(posiciones[i][0]*velocidades[i][1]-posiciones[i][1]*velocidades[i][0]);
    }

    return momento;
}

/*Función integralangulo. Añade el incremento de posición angular a cada planeta. Si el planeta acaba de completar 
la vuelta, cambia el valor del booleano vueltadada*/
void integralangulo(double angulos[], double posiciones[][2], double velocidades[][2], bool vueltadada[]) {

    double radio2;
    double producto;
    
    // Sumo el incremento de posición angular. Para el Sol (i=0) no
    for(int i=1; i<N; i++) {
        if(!vueltadada[i]) {
            radio2=posiciones[i][0]*posiciones[i][0]+posiciones[i][1]*posiciones[i][1];;
            producto=posiciones[i][0]*velocidades[i][1]-velocidades[i][0]*posiciones[i][1];
            angulos[i]+=h*producto/radio2;

            // Si el ángulo supera los 2pi radianes, ha dado la vuelta
            if(abs(angulos[i]>=6.2832)) vueltadada[i]=true;
        }
        
        
    }

    return;

}

/*Función rVerlet. Calcula la nueva posición r a partir de los parámetros r, v y a anteriores. h es el paso temporal*/
double rVerlet(double r, double v, double a) {
    return r+h*v+h*h*a/2;
}


/*Función vVerlet. Calcula la nueva velocidad v a partir de los parámetros v, a anteriores, así como de la aceleración
del mismo paso, a2 h es el paso temporal*/
double vVerlet(double v, double a, double a2) {
    return v+h*(a+a2)/2;
}