//void gammp(double a,double x)

 /******************************************************************** 

  ERRORF 

  this module uses the error function definitions as defined in 
Numerical recipes, Press et al,1992. The error function (erf(y)) is the integral
of the Gaussian distribution and as such yields the discrete surface
areas (surf_factor) over which the total nett sediment has to be distributed.
This is used in the deltaic domain.

NB the errorfunction is used from y=0 to y=2, as an approximation of unity
by definition y-->infinity becomes unity only 
 

**********************************************************************/

double erf(double x)

{
	double gammp(double a,double x);

	return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}

