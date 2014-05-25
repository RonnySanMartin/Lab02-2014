/*
*	DOBLADO ECC CON APOYO DE COORDENADAS JACOBIANAS
*
*	Autor: Ronny A. San Martin Rodriguez - LCC - USACH
*	
*	gcc doubling.c -o doub -lgmp
*	./doub pXXX
*
*	pXXX = Curva eliptica (192, 224, 256, 384, 521) segun estandar NIST
*
*	Cada archivo se compone del n primo y dos puntos (coordenas X e Y c/u)
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "time.h"		//Libreria para calculo del tiempo

#define millon 1000000	//Divide al tiempo gastado para obtenerlo en ms
#define maxchar 160 	//Maximo de caracteres a leer en variable cadena

int main(int argc, char *argv[])
{
	FILE *archivo;			//Archivo fuente con los puntos a sumar
	char cadena[160];		//Cadena auxiliar para pasar de archivo a mpz_t
	time_586 inicio, fin;	//Tiempos
	double tiempo_gastado;

	mpz_t P, b;						//Curva (P primo, b de curva)
	mpz_t x1, y1;					//Punto Entrante (z1=1, calculos triviales)
	mpz_t x3, y3, z3, alpha, beta;	//doblado en coordenadas Jacobianas
	mpz_t x2, y2, z2;				//Punto Resultante
	mpz_t aux1, aux2; 				//Auxiliares

	mpz_inits(P,b, x1,y1, x2,y2,z2, x3,y3,z3, alpha,beta, aux1,aux2,NULL);

	archivo = fopen(argv[1], "r");		//Lectura desde archivo fuente
	if(archivo && argc == 2)
	{
		fgets(cadena, maxchar, archivo);	//Asignacion de variables
		mpz_set_str(P, cadena, 10);
		fgets(cadena, maxchar, archivo);
		mpz_set_str(b, cadena, 10);
		fgets(cadena, maxchar, archivo);
		mpz_set_str(x1, cadena,  10);
		fgets(cadena, maxchar, archivo);
		mpz_set_str(y1, cadena,  10);

		fclose(archivo);

		system("clear");

		time(inicio);//INICIO ALGORITMO

		//Calculo de ALPHA
		mpz_pow_ui(alpha, x1, 2);	//alpha = x1^2
		mpz_sub_ui(alpha, alpha, 1);//alpha = x1^2 - z1^4 (z1^4 = 1)
		mpz_mul_ui(alpha, alpha, 3);//alpha = 3*(x1^2-z1^4)
		mpz_mod(alpha,alpha,P);		//reduccion (modulo)

		//Calculo de BETA
		mpz_mul_ui(beta, x1, 4);	//beta = 4*x1
		mpz_mul(beta, beta, y1);	//beta = 4*x1*y1
		mpz_mul(beta, beta, y1);	//beta = 4*x1*y1*y1
		mpz_mod(beta, beta, P);		//Reduccion

		//Calculo de Z3
		mpz_mul_ui(z3, y1, 2);		//z3 = 2*y1*z1 (z1 = 1, trivial, no se incluye)
		mpz_mod(z3, z3, P);			//Reduccion

		//calculo de X3
		mpz_pow_ui(x3, alpha, 2);	//x3 = alpha^2
		mpz_submul_ui(x3, beta, 2);	//x3 = alpha^2 - 2*beta
		mpz_mod(x3, x3, P);			//Reduccion

		//calculo de Y3
		mpz_pow_ui(aux1, y1, 4);	//aux1 = y1^4
		mpz_sub(y3, beta, x3);		//y3   = beta-x3
		mpz_mul(y3,alpha,y3);		//y3   = alpha*(beta-x3)
		mpz_submul_ui(y3, aux1, 8);	//y3   = alpha*(beta-x3)-8*y1^4
		mpz_mod(y3, y3, P);			//Reduccion

		//Transformacion jacobiano a affine
		mpz_invert(aux1,z3,P);		//aux1 = z3^(-1) mod P
		mpz_mul(z2,aux1,z3);		//z2   = z3*z3^(-1) = 1
		mpz_pow_ui(aux2,aux1,3);	//aux2 = (z3^(-1))^3
		mpz_pow_ui(aux1,aux1,2);	//aux1 = (z3^(-1))^2
		mpz_mul(x2, x3, aux1);		//x2   = x3*(z3^(-1))^2
		mpz_mul(y2, y3, aux2);		//y2   = y3*(z3^(-1))^3
		
		mpz_mod(x2, x2, P);			//Reduccion
		mpz_mod(y2, y2, P);
		mpz_mod(z2, z2, P);

		time(fin);				//FIN ALGORITMO
		tiempo_gastado = time_diff(fin,inicio)/millon;	//Calculo tiempo de ejecucion
		
		//IMPRESION DE DATOS:
		printf("DATOS DE CURVA\n");
		gmp_printf("Primo = %Zd\n", P);
		gmp_printf("b     = %Zd\n", b);
		printf("\nPUNTOS ENTRANTES P y Q\n");
		gmp_printf("x1    = %Zd\n", x1);
		gmp_printf("y1    = %Zd\n", y1);
		gmp_printf("z1    = 1\n");
		printf("\nDATOS DE COORDENADAS JACOBIANAS\n");
		gmp_printf("alpha = %Zd\n",alpha);
		gmp_printf("beta  = %Zd\n",beta);
		gmp_printf("x3    = %Zd\n", x3);
		gmp_printf("y3    = %Zd\n", y3);
		gmp_printf("z3    = %Zd\n", z3);
		printf("\nRESULTADO 2P = (X,Y,Z) (affine)\n");
		gmp_printf("X  = %Zd\n", x2);
		gmp_printf("Y  = %Zd\n", y2);
		gmp_printf("Z  = %Zd\n", z2);

		printf("Tiempo gastado: %.20f\n\n", tiempo_gastado);
	}
	else
	{
		printf("ERROR: uso ./doub pXXX (pXXX = archivo con puntos: p192, p224, p256, p384, p521; segun curva)\n");
	}

	mpz_clears(P,b, x1,y1, x2,y2,z2, x3,y3,z3, alpha,beta, aux1,aux2, NULL);
}