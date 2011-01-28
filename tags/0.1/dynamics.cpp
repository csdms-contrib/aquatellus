// calculates in 2D erosion / sedimentation flux of the longitudinal profile
// which is extracted from the 3D matrix for this specific time step.

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <fstream>
#include <iterator>
#include "main.h"
#include "aquatellus.h"
#include "Space.h"
#include "flowpath.h"

int i = 0;
int k, l, m, n;

double streampower;
double discharge_def = 200;
double discharge_shapef = 0.005;
double waterdepth = 0;

double k_er_fluv = 0.00005;     //0.0000050; //used in experi 1;
double k_er_marine = 0.00001;

// these travel distances reflect the longitudinal travel distance
double traveldist_fluvial[] = { 220, 230, 250, 270 };   // from coarse to fine
double traveldist_marine[] = { 80, 160, 180, 220 };     //from fine to coarse{50,60,80,120};

double traveldist;

double erosion_power;
double sediment_load_def = 0.81;        // the given input sediment load at t=0
double initload;                // the input sediment load at t=0 corrected for climate fluctuation

double settle_rate;
double **sedflux;


flowpath::iterator
determine_coastline_cell (flowpath * current_flowpath, double sealevel,
                          long sim_time)
{
  flowpath::iterator i;
  for (i = current_flowpath->begin (); i != current_flowpath->end (); i++)
  {
    double spaceHeight = get_topo_height (sim_time, i->row, i->col);
    if (spaceHeight < sealevel)
      return i;
  }
  return i;
}

double *
determine_discharge_along_flowpath (flowpath * current_flowpath,
                                    flowpath::iterator coastline_cell_index,
                                    double sealevel, long sim_time)
{
  //DISCHARGE DETERMINATION
  double *discharge;
  if ((discharge =
       (double *) calloc (current_flowpath->size () + 1,
                          sizeof (double))) == NULL)
    printf ("error in initializing discharge array");

  discharge[0] = discharge_def * climate_factor;

  flowpath::iterator i = current_flowpath->begin ();
  i++;                          // discharge[0] is already defined (see above)

  int j = 1;
  // onshore
  for (; i != coastline_cell_index && i != current_flowpath->end (); i++, j++)
  {
    discharge[j] = discharge[j - 1];
  }


  for (; i != current_flowpath->end (); i++, j++)       // offshore
  {
    //erosion decreases with increasing waterdepth
    //this is mimicked by manipulating the discharge
    double waterdepth = sealevel - get_topo_height (sim_time, i->row, i->col);


    discharge[j] =
      discharge[j - 1] -
      (discharge_shapef * dx * waterdepth * discharge[j - 1]);

    if (discharge[j] < 0.0)
    {
      discharge[j] = 0.0;
    }

  }

  return discharge;
}


double *
determine_erosion_rate_along_flowpath (flowpath * current_flowpath,
                                       flowpath::
                                       iterator coastline_cell_index,
                                       double *discharge, long sim_time)
{

  // EROSION DETERMINATION (GRAINSIZE INDEPENDENT)
  double erosion_power = 1.0;
  double slope;

  double *dh_erosion;
  if ((dh_erosion =
       (double *) calloc (current_flowpath->size () + 1,
                          sizeof (double))) == NULL)
    printf ("error in initializing dh_erosion array");

  int j = 0;
  // this flowpath iterator is necessary to look one cell forward
  // this is necessary to calculate the local slope
  flowpath::iterator nextCell = current_flowpath->begin ();
  if (nextCell != current_flowpath->end ())
    nextCell++;

  for (flowpath::iterator i = current_flowpath->begin ();
       i != current_flowpath->end (); i++, j++, nextCell++)

  {
    //erosion is slope independent in the marine domain
    // dit effect heb je nu natuurlijk niet!!
    //if(i == coastline_cell_index) 
    //{ 
    // dit was erosion_power = 0.0; 
    //maar dat heeft niet zoveel zin als slope een getal kleiner dan 0 is
    // dan is het netto effect alleen maar dat de erosion depth kleiner wordt in
    // het fluviatiele gedeelte tov het mariene gedeelte
    //      erosion_power = 1.0;
    //}


    if (nextCell != current_flowpath->end ())
    {
      slope =
        (get_topo_height (sim_time, i->row, i->col) -
         get_topo_height (sim_time, nextCell->row, nextCell->col)) / dx;
    }
    // in this way the lower boundary condition for slope(at the end of the grid) is set 
    // last slope = previous slope


    //erosion depends with a power function on local slope and discharge
    double streampower = discharge[j] * pow (slope, erosion_power);

    // the erosion constants in the equation are different
    // for the fluvial and marine environment
    double erosion_rate;
    if (i != coastline_cell_index)
    {
      erosion_rate = k_er_fluv * streampower;

    }
    else
    {
      erosion_rate = k_er_marine * streampower;

    }

    dh_erosion[j] = dt * erosion_rate;


  }
  return dh_erosion;
}




double **
dig_erosion_from_space (flowpath * current_flowpath, double *dh_erosion,
                        long sim_time)
{
  double **sedflux_ero;
  if ((sedflux_ero =
       (double **) calloc (current_flowpath->size () + 1,
                           sizeof (double *))) == NULL)
    printf ("error in initializing erosional sedflux array");

  for (int k = 0; k < current_flowpath->size () + 1; k++)
  {
    if ((sedflux_ero[k] =
         (double *) calloc (number_of_gscl, sizeof (double))) == NULL)
      printf ("error in initializing erosional sedflux array");
  }

  //erosion is substracted from the space array

  int j = 0;
  for (flowpath::iterator i = current_flowpath->begin ();
       i != current_flowpath->end (); i++, j++)
  {
    bool done = false;
    for (long m = sim_time - 1; m >= 0 && !done; m--)
    {
      int height = space[m][i->row][i->col][0];
      if (height != 0)
      {
        if (height - dh_erosion[j] < 0)
        {
          // Remove all sediment from space for this specific time
          for (int l = 1; l <= number_of_gscl; l++)
          {
            // is het hier wel *(dx/dt)?
            sedflux_ero[j][l - 1] += space[m][i->row][i->col][l] * (dx / dt);
            space[m][i->row][i->col][l] = 0;
          }
          dh_erosion[j] -= space[m][i->row][i->col][0];
          space[m][i->row][i->col][0] = 0;
        }
        else
        {
          // Remove dh_erosion from stored amount space
          for (int l = 1; l <= number_of_gscl; l++)
          {
            double grain_percentage =
              grain_percent (m, i->row, i->col, l) * dh_erosion[j];
            sedflux_ero[j][l - 1] += grain_percentage * (dx / dt);
            space[m][i->row][i->col][l] -= grain_percentage;
          }
          space[m][i->row][i->col][0] -= dh_erosion[j];
          done = true;
        }
      }
    }
  }
  return sedflux_ero;
}

double **
sedimentation_rate_along_flowpath (flowpath * current_flowpath,
                                   flowpath::iterator coastline_cell_index,
                                   double **sedflux_ero)
{

  double **sedflux_depo;
  if ((sedflux_depo =
       (double **) calloc (current_flowpath->size () + 1,
                           sizeof (double *))) == NULL)
    printf ("error in initializing depositional sedflux array");

  for (int k = 0; k < current_flowpath->size () + 1; k++)
  {
    if ((sedflux_depo[k] =
         (double *) calloc (number_of_gscl, sizeof (double))) == NULL)
      printf ("error in initializing depositional sedflux array");
  }

// SEDIMENTATION DETERMINATION (GRAINSIZE SPECIFIC)

  // INITIAL SEDIMENT LOAD DETERMINATION
  initload = sediment_load_def * climate_factor;

  for (n = 0; n < number_of_gscl; n++)
  {
    // calculate the distribution of the total load 
    // over the different grainsize classes
    sedflux_depo[0][n] = sedflux_ero[0][n] + (initload * sed_cont_pct[n]);
  }


  flowpath::iterator i = current_flowpath->begin ();
  i++;                          // sedflux_depo[0] is already defined (see above)
  int j = 1;
  double settle_rate;
  // onshore
  for (; i != coastline_cell_index && i != current_flowpath->end (); i++, j++)
  {
    for (n = 0; n < number_of_gscl; n++)
    {
      traveldist = traveldist_fluvial[n];
      settle_rate =
        (sedflux_depo[j - 1][n] + sedflux_ero[j - 1][n]) / traveldist;
      sedflux_depo[j][n] = sedflux_depo[j - 1][n] - (dx * settle_rate);
      //cout<<" j "<<j<<" sedfluxdepo = "<<sedflux_depo[j][n]<<endl;
    }

  }

  for (; i != current_flowpath->end (); i++, j++)       // offshore              
  {

    for (n = 0; n < number_of_gscl; n++)
    {
      traveldist = traveldist_marine[n];
      settle_rate =
        (sedflux_depo[j - 1][n] + sedflux_ero[j - 1][n]) / traveldist;
      sedflux_depo[j][n] = sedflux_depo[j - 1][n] - (dx * settle_rate);
    }
  }

  return sedflux_depo;
}
