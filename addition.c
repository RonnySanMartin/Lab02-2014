/*
*	ADICION EN ECC CON APOYO DE COORDENADAS JACOBIANAS
*
*	Autor: Ronny A. San Martin Rodriguez - LCC - USACH
*	
*	gcc addition.c -o add -lgmp
*	./add
*
* 	Curva: y^2 = x^3+ax+b
* 	P = 6277101735386680763835789423207666416083908700390324961279
* 
*	a = -3
* 	b = 2455155546008943817740293915197451784769108058161191238065
*
* 	Punto P = (X1,Y1,Z1):
* 	X1 = 602046282375688656758213480587526111916698976636884684818
* 	Y1 = 174050332293622031404857552280219410364023488927386650641
* 	Z1 = 1
*
*	Punto Q = (X2,Y2,Z2):
*	X2 = 842308855964553312804811126459961974190944517868650739565
*	Y2 = 4419613963320162418875315725715541668384579137277197535329
*	Z2 = 1
*
*	P+Q = (X4,Y4,Z4):
*	X4 = 2677124060146511202177285576496406332332084632597795659103
*	Y4 = 108464341913100350370050211857368086435213456474211256786
*	Z4 = 1
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "time.h"		//Libreria para calculo del tiempo

#define millon 1000000	//Divide al tiempo gastado para obtenerlo en ms

int main(void)
{
	time_586 inicio, fin;
	double tiempo_gastado;

	mpz_t P, b;						//Curva (P primo, b de curva)
	mpz_t x1, y1, x2, y2;			//Punto Entrante (z1=1, calculos triviales)
	mpz_t x3, y3, z3, alpha, beta;	//Adicion en coordenadas Jacobianas
	mpz_t x4, y4, z4;				//Punto Resultante
	mpz_t aux1, aux2; 				//Auxiliares

	mpz_inits(P,b, x1,y1, x2,y2, x3,y3,z3, x4,y4,z4, alpha,beta, aux1,aux2,NULL);

	mpz_set_str(P,  "6277101735386680763835789423207666416083908700390324961279", 10);
	
	mpz_set_str(x1, "602046282375688656758213480587526111916698976636884684818",  10);
	mpz_set_str(y1, "174050332293622031404857552280219410364023488927386650641",  10);
	
	mpz_set_str(x2, "842308855964553312804811126459961974190944517868650739565",  10);
	mpz_set_str(y2, "4419613963320162418875315725715541668384579137277197535329", 10);	

	system("clear");

	time(inicio);//INICIO ALGORITMO

	//Calculo de ALPHA
	mpz_sub(alpha, y2, y1);		//alpha = z1^3*y2 - z2^3*y1 (z1 = z2 = 1, se omite)
	mpz_mod(alpha, alpha, P);	//Reduccion

	//Calculo de Beta
	mpz_sub(beta, x2, x1);		//beta = z1^2*x2 - z2^2*x1 (z1 = z2 = 1, se omite)
	mpz_mod(beta, beta, P);		//Reduccion

	//Calculo de z3
	mpz_set(z3, beta);		//z3 = z1*z2*beta (z1 = z2 = 1, se omite)
							//no se reduce ya que beta se encuentra reducido

	//Calculo de x3
	mpz_pow_ui(x3, alpha, 2);	//x3   = alpha^2
	mpz_pow_ui(aux1, beta, 2);	//aux1 = beta^2
	mpz_submul(x3, aux1, beta);	//x3   = alpha^2 - beta^3
	mpz_mul(aux1, aux1,x1);		//aux1 = z2^2*x1*beta^2 (z2 = 1, se omite)
	mpz_submul_ui(x3, aux1, 2);	//x3   = alpha^2 - beta^3 -2*z2^2*x1*beta^2 (z2 = 1, se omite)
	mpz_mod(x3, x3, P);			//Reduccion

	//calculo de y3
	mpz_sub(y3, aux1, x3);		//y3   = z2^2*x1*beta^2 - x3
	mpz_mul(y3, alpha, y3);		//y3   = aplha*(z2^2*x1*beta^2 - x3)
	mpz_pow_ui(aux1, beta, 3);	//aux1 = z2^3*beta^3
	mpz_submul(y3, aux1, y1);	//y3   = aplha*(z2^2*x1*beta^2 - x3) - z2^2*y1*beta^3
	mpz_mod(y3, y3, P);			//Reduccion

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
	printf("DATOS DE COORDENADAS JACOBIANAS\n\n");
	gmp_printf("alpha = %Zd\n", alpha);
	gmp_printf("beta  = %Zd\n", beta);
	gmp_printf("x3    = %Zd\n", x3);
	gmp_printf("y3    = %Zd\n", y3);
	gmp_printf("z3    = %Zd\n", z3);
	printf("\nRESULTADO P+Q = (X,Y,Z) (affine)\n\n");
	gmp_printf("X  = %Zd\n", x4);
	gmp_printf("Y  = %Zd\n", y4);
	gmp_printf("Z  = %Zd\n", z4);

	printf("Tiempo gastado: %.20f\n\n", tiempo_gastado);

	mpz_clears(P,b, x1,y1, x2,y2, x3,y3,z3, x4,y4,z4, alpha,beta, aux1,aux2, NULL);
}