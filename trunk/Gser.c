#include <stdio.h>
#include <math.h>

#define ITMAX 100
#define EPS 3.0e-7
double gammln(double xx);

void gser(double * gamser,double  a,double  x,double * gln)
{
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) printf("x less than 0 in routine GSER\n");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		printf("a too large, ITMAX too small in routine GSER\n");
		return;
	}
}

#undef ITMAX
#undef EPS
