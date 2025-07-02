#include <R.h>
 #include <math.h>
 void g_Notal_deriv_ojp2oqzu ( double * x, double * y, double * p, int * n, int * k, int * l ) {
 for(int i = 0; i< *n; i++) {
 y[0+i**l] = (1.0)*(x[4+i**k]) ;
y[1+i**l] = (1.0)*(x[7+i**k]) ;
y[2+i**l] = (1.0)*(x[10+i**k]) ;
y[3+i**l] = (1.0)*(x[13+i**k]) ;
y[4+i**l] = (1.0)*(x[16+i**k]) ;
y[5+i**l] = (1.0)*(x[19+i**k]) ;
y[6+i**l] = (1.0)*(x[22+i**k]) ;
y[7+i**l] = (1.0)*(x[25+i**k]) ;
y[8+i**l] = (1.0)*(x[28+i**k]) ; 
}
}