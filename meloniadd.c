/*
*	ADICION EN ECC CON APOYO DE FORMULAS DE MELONI
*
*	Autor: Ronny A. San Martin Rodriguez - LCC - USACH
*	
*	gcc meloniadd.c -o madd -lgmp
*	./madd pXXX
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

	mpz_t P, b;				//Curva (P primo, b de curva)
	mpz_t x1, y1, x2, y2; 	//Puntos entrantes (z=1, se omite)
	mpz_t mA, mB, mC, mD;	//Variables para F. Meloni
	mpz_t x3, y3, z3;		//Punto resultante
	mpz_t x4, y4, z4;
	mpz_t aux1, aux2;		//Auxiliares

	mpz_inits(P,b, x1,y1, x2,y2, x3,y3,z3, x4,y4,z4, mA,mB,mC,mD, aux1,aux2, NULL);

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
		fgets(cadena, maxchar, archivo);
		mpz_set_str(x2, cadena,  10);
		fgets(cadena, maxchar, archivo);
		mpz_set_str(y2, cadena, 10);	

		fclose(archivo);

		system("clear");

		time(inicio);//INICIO ALGORITMO

		//Calculo de mA
		mpz_sub(mA, x2, x1);
		mpz_pow_ui(mA, mA, 2);
		mpz_mod(mA, mA, P);		//Reduccion

		//Calculo de mB
		mpz_mul(mB, x1, mA);
		mpz_mod(mB, mB, P);		//Reduccion

		//Calculo de mC
		mpz_mul(mC, x2, mA);
		mpz_mod(mC, mC, P);		//Reduccion

		//Calculo de mD
		mpz_sub(mD, y2, y1);
		mpz_pow_ui(mD, mD, 2);
		mpz_mod(mD, mD, P);		//Reduccion

		//Calculo de X3
		mpz_sub(x3, mD, mB);
		mpz_sub(x3, x3, mC);
		mpz_mod(x3, x3, P);		//Reduccion

		//Calculo de Y3
		mpz_sub(aux1, y2, y1);
		mpz_sub(y3, mB, x3);
		mpz_mul(y3, y3, aux1);
		mpz_sub(aux1, mC, mB);
		mpz_submul(y3, aux1, y1);
		mpz_mod(y3, y3, P);		//Reduccion

		//Calculo de Z3
		mpz_sub(z3, x2, x1);
		mpz_mod(z3, z3, P);		//Reduccion

		//Transformacion jacobiano a affine
		mpz_invert(aux1,z3,P);		//aux1 = z3^(-1) mod P
		mpz_mul(z4,aux1,z3);		//z2   = z3*z3^(-1) = 1
		mpz_pow_ui(aux2,aux1,3);	//aux2 = (z3^(-1))^3
		mpz_pow_ui(aux1,aux1,2);	//aux1 = (z3^(-1))^2
		mpz_mul(x4, x3, aux1);		//x2   = x3*(z3^(-1))^2
		mpz_mul(y4, y3, aux2);		//y2   = y3*(z3^(-1))^3
		
		mpz_mod(x4, x4, P);			//Reduccion
		mpz_mod(y4, y4, P);
		mpz_mod(z4, z4, P);


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
		gmp_printf("x2    = %Zd\n", x2);
		gmp_printf("y2    = %Zd\n", y2);
		gmp_printf("z2    = 1\n");
		printf("\nDATOS DE FORMULA DE MELONI\n");
		gmp_printf("A     = %Zd\n", mA);
		gmp_printf("B     = %Zd\n", mB);
		gmp_printf("C     = %Zd\n", mC);
		gmp_printf("D     = %Zd\n", mD);
		printf("\nRESULTADO P+Q = (X',Y',Z')\n");
		gmp_printf("X'    = %Zd\n", x4);
		gmp_printf("Y'    = %Zd\n", y4);
		gmp_printf("Z'    = %Zd\n", z4);
		printf("\nTiempo gastado: %.20f\n\n", tiempo_gastado);
	}
	else
	{
		printf("ERROR: uso ./madd pXXX (pXXX = archivo con puntos: p192, p224, p256, p384, p521; segun curva)\n");
	}

	mpz_clears(P,b, x1,y1, x2,y2, x3,y3,z3, x4,y4,z4, mA,mB,mC,mD, aux1,aux2, NULL);

}