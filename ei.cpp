/******************************************************

  EXPINT this function is based on numerical recipes
  evaluates the exponential integral Ei(x)

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
	int i,nval;
	//float val, x;
	float x;
	FILE *fp, *out;

	double ei(float x);
	out = fopen("ei_out.dat", "w");

	if ((fp=fopen( "ei_inp.dat", "r")) ==NULL)
		printf("file nor found\n");
	fgets(txt, MAXSTR,fp);
	while (strncmp(txt, "Exponential Integral Ei", 23)) {
		fgets(txt, MAXSTR,fp);
		if(feof(fp)) printf("data not found in file\n");
		}

	fscanf(fp, "%d %*s", &nval);
	printf("\n%s\n",txt);
	//printf("%5s %12s %11s \n","x","actual","ei(x)");
	printf("%5s %11s \n","x","ei(x)");
	for (i=1; i<=nval; i++) {
	//	fscanf(fp," %f %f",&x,&val);
		fscanf(fp," %f",&x);
		//printf("%6.2f %12.6f %12.6f\n",x,val,ei(x));
		//fprintf(out,"%6.2f %12.6f %12.6f\n",x,val,ei(x));
		printf("%6.2f %12.6f\n",x, ei(x));
		fprintf(out,"%6.2f %12.6f\n",x,ei(x));
	}
	fclose(fp);
	fclose(out);
	return 0;
}
*/

double ei(float x)
{
	int k;
	double fact,prev,sum,term;

	
	if (x<=0.0)	printf("bad arguments in ei\n");
	if (x<FPMIN)return log(x) + EULER;
	if (x<= -log(EPS)) {
		sum = 0.0;
		fact = 1.0;
		for (k=1; k<=MAXIT; k++) {
			fact *= x/k;
			term = fact/k;
			sum+= term;
			if (term < EPS*sum) break;

		}
		 
		if (k> MAXIT) printf("series failed in ei\n");
		return sum +log(x)+EULER;
	}
	else{
		sum =0.0;
		term =1.0;
		for (k=1; k<=MAXIT; k++) {
			prev =term;
			term *= k/x;
			if (term < EPS) break;
			// since the final sum is greater than one
			// term itself approximates the relative error

			if(term< prev) sum+= term;
			else{
				sum-= prev;
				break;
			}
		}
		return exp(x) * (1.0 + sum)/x;
	}
}




