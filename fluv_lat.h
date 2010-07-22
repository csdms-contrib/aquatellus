#ifndef FLUV_LAT_H_
#define FLUV_LAT_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <iostream>


#include "main.h"
#include "Space.h"
#include "flowpath.h"
#include "inputcontrols.h"
#include "dynamics.h"
#include "output.h"
#include "ei.h"



#define MAXIT 100               // maximum number of iterations
#define EULER 0.5772156649      //Euler's constant gamma
#define FPMIN 1.0E-30           // close to smallest representable float
#define EPS1 1.0E-07            //desired relative error
#define MAXSTR 80               // for decapping the file headers

typedef struct
{
  int left;
  int right;
} COORPAIR;


double *calculate_firstder (flowpath * current_flowpath);
double *spline (flowpath * current_flowpath, double *firstder);
double *calculate_curvature (flowpath * current_flowpath, double *firstder,
                             double *y2);
double *calculate_asymmetry (flowpath * current_flowpath, double *curvature);
COORPAIR storage (flowpath::iterator i, double *discharge, int j,
                  long sim_time);
double expint (int n, double x);
void fluv_lat (flowpath * current_flowpath,
               flowpath::iterator coastline_cell_index, long sim_time,
               double *discharge, double **sedflux_depo);


#endif
