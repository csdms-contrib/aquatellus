/******************************************************

  EXPINT this function is based on numerical recipes
  evaluates the exonential integral En(x)

  ******************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXIT 100			// maximum number of iterations
#define EULER 0.5772156649	//Euler's constant gamma
#define FPMIN 1.0E-30		// close to smallest representable float
#define EPS 1.0E-07			//desired relative error
#define MAXSTR 80			// for decapping the file headers

/*
int main (void)

{
	char txt[MAXSTR];
	int i,nval,n;
	float val, x;
	FILE *fp, *out;

	float expint(int n,float x);
	out = fopen("expint_out.dat", "w");

	if ((fp=fopen( "expdat2.txt", "r")) ==NULL)
		printf("file not found\n");
	fgets(txt, MAXSTR,fp);
	while (strncmp(txt, "Exponential Integral En", 23)) {
		fgets(txt, MAXSTR,fp);
		if(feof(fp)) printf("data not found in file\n");
		}

	fscanf(fp, "%d %*s", &nval);
	printf("\n%s\n",txt);
	printf("%4s %20s %21s \n", "n","x","expint(n,x)");
	for (i=1; i<=nval; i++) {
		fscanf(fp,"%d %f ", &n,&x);
		printf("%4d %8.2f %18.6e\n", n,x,expint(n,x));
		fprintf(out,"%4d %8.2f %18.6e\n", n,x,expint(n,x));
	}
	fclose(fp);
	fclose(out);
	return 0;
}
*/
/*
double expint(int n,double x)
{
	int i,ii,nm1;
	double a,b,c,d,del,fact,h,psi,ans;

	nm1 = n-1;
	if (n<0||x<0.0||(x==0.0&&(n==0||n==1)))
		printf("bad arguments in expint\n");
	else {
		if (n==0) ans=exp(-x)/x;	//special case
		else {
			if(x==0.0) ans =1.0/nm1; //another special case

			else {

				if (x>1.0) {		//Lentz's algorithm
					b = x+n;
					c = 1.0/FPMIN;
					d = 1.0/b;
					h = d;
					for(i=1; i< MAXIT; i++) {
						a = -i * (nm1+i);
						b += 2.0;
						d = 1.0/(a*d+b);
						c = b+a/c;
						del = c*d;
						h*=del;
						if (fabs(del-1.0) < EPS){
							ans = h*exp(-x);
							return ans;
						}
					}
					printf("continued fraction failed in expint\n");

				} else {
					ans = (nm1!= 0?1.0/nm1 : -log(x)-EULER);
					fact = 1.0;
					for (i=1; i <MAXIT; i++) {
						fact *= -x/i;
						if(i!=nm1) del = -fact/(i-nm1);
						else {
							psi = -EULER;
							for (ii =1; ii<=nm1; ii++) psi += 1.0/ii;
							del = fact * (-log(x)+ psi);
						}
						ans +=del;
						if (fabs(del) < fabs(ans) *EPS) return ans;

					}
					printf("series failed in expint\n");
				}
			}
		}
	}
	return ans;
}


*/