/*
 *  aquatellusreadinput.cpp
 *
 *
 *  Created by Irina Overeem, September 2011.
 *  Copyright 2011 Univ of CO. All rights reserved.
 *
 */

#include "main.h"
#include "aquatellus.h"
#include "aquatellusreadinput.h"
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>

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
double sea_lvl_ref = 80;        // average level around which fluctuations occur, reference sea level at t=0
double sea_lvl_amp = 0;         // amplitude in m
double sea_lvl_prd = 100;       // period of the sine function in years

void
aquatellusreadinput ()
{
  FILE *inp;

  char chs[150];

  int jj = 0;

  if ((inp = fopen ("aquatellusinput.rtf", "r")) == NULL)
  {
    printf ("Error in opening file\n");
    fflush (stdout);
  }
  else
  {
    printf ("file successfully opened\n");
    fflush (stdout);

    // read header
    fgets (chs, 150, inp);
    fgets (chs, 150, inp);
    fgets (chs, 150, inp);

    //SIMULATION TIME SETTINGS
    // read simulation time
    fscanf (inp, "%ld", &end_of_times);
    printf ("simulation time = %ld\n", end_of_times);
    fflush (stdout);
    fgets (chs, 150, inp);

    // read time step
    fscanf (inp, "%d", &dt);
    printf ("timestep = %d\n", dt);
    fflush (stdout);
    fgets (chs, 150, inp);

    // READ GRID DIMENSIONS PARAMETERS

    // read number of grid rows
    fscanf (inp, "%d", &number_of_rows);
    printf ("number  of grid rows = %d\n", number_of_rows);
    fflush (stdout);
    fgets (chs, 150, inp);

    // read number of grid columns
    fscanf (inp, "%d", &number_of_colums);
    printf ("number  of grid colums = %d\n", number_of_colums);
    fflush (stdout);
    fgets (chs, 150, inp);

    //read dx and dy for grid
    fscanf (inp, "%d", &dx);
    printf ("dx = %d\n", dx);
    fflush (stdout);
    fgets (chs, 150, inp);
    fscanf (inp, "%d", &dy);
    printf ("dy = %d\n", dy);
    fflush (stdout);
    fgets (chs, 150, inp);

    // READ GRID GEOMETRY

    // read starting height at upstream end of floodplain
    fscanf (inp, "%lf", &initial_height);
    printf ("height at upstream end of floodplain = %lf\n", initial_height);
    fflush (stdout);
    fgets (chs, 150, inp);

    // read initial gradient in river floodplain
    fscanf (inp, "%lf", &initial_gradient);
    printf ("river floodplain gradient = %lf\n", initial_gradient);
    fflush (stdout);
    fgets (chs, 150, inp);

    // read initial gradient for bathymetry
    fscanf (inp, "%lf", &offshore_gradient);
    printf ("offshore gradient = %lf\n", offshore_gradient);
    fflush (stdout);
    fgets (chs, 150, inp);

    // READ SEDIMENT PARAMETERS

    // read number of grain size classes
    fscanf (inp, "%d", &number_of_gscl);
    printf ("number  of grain size classes = %d\n", number_of_gscl);
    fflush (stdout);
    fgets (chs, 150, inp);

    // read grain size class in mm
    for (jj = 0; jj < number_of_gscl; jj++)
    {
      fscanf (inp, "%lf", &grain_size[jj]);
      printf ("grainsize %d = %lf\n", jj + 1, grain_size[jj]);
      fflush (stdout);
    }
    fgets (chs, 150, inp);

    // read percentage of sediment in each grain size class
    for (jj = 0; jj < number_of_gscl; jj++)
    {
      fscanf (inp, "%lf", &sed_cont_pct[jj]);
      printf ("percentage of grainsize %d = %lf\n", jj + 1, sed_cont_pct[jj]);
      fflush (stdout);
    }
    fgets (chs, 150, inp);

    //read travel distance in fluvial transport for each grain size class
    for (jj = 0; jj < number_of_gscl; jj++)
    {
      fscanf (inp, "%lf", &traveldist_fluvial[jj]);
      printf ("travel distance in fluvial domain %d = %lf\n", jj + 1,
              traveldist_fluvial[jj]);
      fflush (stdout);
    }
    fgets (chs, 150, inp);

    //read travel distance for marine transport for each grain size class
    for (jj = 0; jj < number_of_gscl; jj++)
    {
      fscanf (inp, "%lf", &traveldist_marine[jj]);
      printf ("travel distance in marine domain %d = %lf\n", jj + 1,
              traveldist_marine[jj]);
      fflush (stdout);
    }
    fgets (chs, 150, inp);

    // read k_erosion in fluvial domain
    fscanf (inp, "%lf", &k_er_fluv);
    printf ("fluvial erosion proportionality constant = %lf\n", k_er_fluv);
    fflush (stdout);
    fgets (chs, 150, inp);

    // read k_erosion in marine domain
    fscanf (inp, "%lf", &k_er_marine);
    printf ("marine erosion proportionality constant = %lf\n", k_er_marine);
    fflush (stdout);
    fgets (chs, 150, inp);

    // READ CONTROLLING RIVER PARAMETERS
    fscanf (inp, "%lf", &discharge_def);
    printf ("incoming river discharge = %lf\n", discharge_def);
    fflush (stdout);
    fgets (chs, 150, inp);

    fscanf (inp, "%lf", &sediment_load_def);
    printf ("incoming river sediment load = %lf\n", sediment_load_def);
    fflush (stdout);
    fgets (chs, 150, inp);

    // READ CONTROLLING SEA LEVEL PARAMETERS

    fscanf (inp, "%lf", &sea_lvl_ref);
    printf ("sea level at t=0 = %lf\n", sea_lvl_ref);
    fflush (stdout);
    fgets (chs, 150, inp);

    fscanf (inp, "%lf", &sea_lvl_amp);
    printf ("sea level amplitude of sine function = %lf\n", sea_lvl_amp);
    fflush (stdout);
    fgets (chs, 150, inp);

    fscanf (inp, "%lf", &sea_lvl_prd);
    printf ("sea level period of sine function = %lf\n", sea_lvl_prd);
    fflush (stdout);
    fgets (chs, 150, inp);

    // close the input file
    fclose (inp);

  }

}
