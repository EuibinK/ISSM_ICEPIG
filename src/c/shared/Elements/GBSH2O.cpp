/* \file GBSH2O.cpp
 * brief figure out B of H2O ice for a certain temperature
 *		INPUT function B=GBSH2O(temperature)
 *    	where rigidigty (in s^(1/n)Pa) is the flow law paramter in the flow law sigma=B*e(1/n) (Nye, p2000). 
 */

#include "../io/io.h"
#include <math.h>
#include "../Numerics/types.h"

IssmDouble GBSH2O(IssmDouble temperature){

	/*Coefficients*/
	const IssmPDouble Rg      = 8.3144598; /* J mol^-1 K^-1 */
	const IssmPDouble n       = 1.8;        /*Glen's exponent*/
	
	const IssmPDouble rd      = 1.5e-6;
	const IssmPDouble f       = 0.25;
	const IssmPDouble phi     = 0.01;
	const IssmPDouble A0      = 3.9e-3;      /*s^-1 MPa       */
	const IssmPDouble p       = -1.4;
	IssmDouble d = (13.4e-6/pow(f,0.5))*(rd/1.5e-6)*pow((0.05/phi),0.5);
	IssmDouble A1 = A0*0.5*pow(3,(n+1)/2);
	IssmDouble A2 = A1*exp(-2*n*phi)*pow(d,p);
	
	const IssmPDouble Q       = 49000.;    /*J mol^-1       */

	/*Arrhenius Law*/
	IssmDouble A = A2*exp(-Q/(temperature*Rg));  /*s^-1 MPa   */
	IssmDouble B = 1e6*pow(A,-1/n);                    /*s^(1/n) Pa */

	/*Beyond-melting-point case*/
	if(temperature>=273.15) _printf0_("H2O ICE - GUARANTEED MELTING. Some temperature values are beyond 273.15K.\n");

	/*Return output*/
	return B;
}
