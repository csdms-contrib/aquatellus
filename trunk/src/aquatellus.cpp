///////////////////////////////////////////////////////////////////////////////
// Aqua-Tell-Us - AQUA = WATER, TELLUS = EARTH -
//
// simulates fluvio-deltaic processes in 3D on a geological scale
// erosion and depostion for several grainsize classes is calculated
// the model is focused on system response over long temporal and spatial scales
//
// Author: Irina Overeem, CSDMS, Institute of Arctic and Alpine Research
// University of Colorado, Boulder, Colorado
// irina.overeem@colorado.edu
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "main.h"
#include "aquatellus.h"
#include "Space.h"
#include "flowpath.h"
#include "inputcontrols.h"
#include "dynamics.h"
#include "fluv_lat.h"
#include "del_lat.h"
#include "output.h"
#include "ei.h"
#include "aquatellus_init_func.h"
#include "aquatellusreadinput.h"

#include "bmi.h"

using namespace std;

// time & grid dimensions
/*
long end_of_times = 200;

int dt = 1;

int number_of_time_steps = end_of_times / dt;

int number_of_rows = 1;

int number_of_colums = 1;

int dx = 1;

int dy = 1;

double initial_height = 0.0;    // starting height at the upstream end of the grid

double initial_gradient = 0.0;  // initial gradient in longitudinal direction in river floodplain

double offshore_gradient = 0.0; // initial gradient in longitudinal direction in offshore floodplain

int number_of_gscl = 3;         // number of grain size classes
double grain_size[] = { 0.0, 0.0, 0.0 };        // grainsizes in mm, coarse to fine
double sed_cont_pct[] = { 0.33, 0.33, 0.33 };   //percentage of each grainsize
double traveldist_fluvial[] = { 220, 230, 270 };        //{ 220, 250, 270, 320 };   // from coarse to fine, set originally { 220, 230, 2 50, 270 }
double traveldist_marine[] = { 50, 60, 120 };   //from fine to coarse {50,60,80,120};

double k_er_fluv = 0.0005;      //0.0000050; //used in experi 1;

double k_er_marine = 0.00001;

double discharge_def = 200;

double sediment_load_def = 0.02;        // the given input sediment load at t=0

///////////////////////////Sea level controls -sine functions//////////////////////////////////
double sea_lvl_ref = 80;        // average level around which fluctuations occur, reference sea level at t=0

double sea_lvl_amp = 0;         // amplitude in m

double sea_lvl_prd = 100;       // period of the sine function in years
*/
// some specifications for output
//int Xpos = 8;                   //position at which Xsection for output is made

int num_outp = 20;              // number of times that topography file gets written away

// X-Y position of simulated wells
int wellpos_row = 45;

int wellpos_col = 85;

double *discharge;

double *dh_erosion;

int **flowpath_array;

double *firstder;

COORPAIR curve_width;

//long sim_time;

//double sealevel;

//double climate_factor;

///////////////////// THE ACTUAL SIMULATION//////////////////////////////////
int
main ()
{

  // INITIALIZATION
  printf ("Reading Aquatellus Input file\n");
  printf ("Initializing DEM and Sediment Grid\n");

  BMI_Initialize (NULL);

  // TIME LOOP

  sim_time = 0;
  while (sim_time < end_of_times)

  {

    //run_aquatellus ();
    BMI_Run_model (NULL);

    if (sim_time >= end_of_times)
    {
      printf ("end of programme\n");
    }
    //printf("simtime=%d\n",sim_time);
    sim_time += dt;

  }     //end of simulation time loop

  //aquatellus_finalize ();
  BMI_Finalize (NULL);

  return 0;

}
