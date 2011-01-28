//////////////////////////////////////////////////////////////////////////
//
//  DEL_LAT this function determines the lateral distribution of sediment
//  in the deltaic domain
//
//  1 calculation of hydrodynamical conditions
//  2 decide bedfriction-plume or turbulent jet plume
//  3 determine lateral width of curve + normalize
//  4 determine lateral sedimentation distribution
//////////////////////////////////////////////////////////////////////////////

#include "Nr.h"
//#include "nrutil.h"


#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <fstream>

#include "main.h"
#include "aquatellus.h"
#include "Space.h"
#include "flowpath.h"
#include "inputcontrols.h"
#include "dynamics.h"
#include"del_lat.h"


void
del_lat (flowpath * current_flowpath, flowpath::iterator coastline_cell_index,
         double *discharge, long sim_time, double **sedflux_depo)
{
#define VISCO  0.8
  int lat_width = 0;
  double b, h, u, b_h_ratio, Re;
  double tan_alpha;

  //dynamically declared for each timestep
  int plume_length = number_of_rows - coastline_cell_index->row;

  double *l_scale;
  l_scale = (double *) calloc (number_of_colums, sizeof (double));

  // to divide the depositional sediment flux over a left and righthand component
  double *sedflux_r;
  sedflux_r = (double *) calloc (number_of_gscl, sizeof (double));

  double *sedflux_l;
  sedflux_l = (double *) calloc (number_of_gscl, sizeof (double));


  int flx = coastline_cell_index->row;  // the row position of the coastline cell in the matrix
  int fly = coastline_cell_index->col;  // the column position of the coastline cell in the matrix

  int j = 0;
  flowpath::iterator i = current_flowpath->begin ();
  i++;
  j++;

  // onshore nothing happens yet, just passing by
  for (; i != coastline_cell_index; i++, j++)
  {
    printf (".");
  }

  if (i == coastline_cell_index && i != current_flowpath->end ())       // at the coastline cell
  {
    printf ("+");
    // this step calculates the hydrodynamical dimensions at the river mouth
    // based on discharge,the formula's are modified from Leopold and Maddock (1953)

    b = pow (discharge[j], 0.5);        //riverwidth
    h = pow (discharge[j], 0.4);        // riverdepth
    u = pow (discharge[j], 0.1);        // river velocity
    b_h_ratio = b / h;          // width/depth ratio
    Re = u * (pow ((h * (b / 2)), 0.5)) / VISCO;        //Reynolds number

    double waterdepth =
      get_topo_height (sim_time, coastline_cell_index->row,
                       coastline_cell_index->col) - sealevel;

    // here the plume type is decided

    if ((Re > 3000) && (h < waterdepth))        //bed friction
    {
      tan_alpha = 0.364;
      printf (" the plume is bedfriction dominated\n");
    }

    else                        //turbulent jet
    {
      tan_alpha = 0.249;
      printf ("the plume is a turbulent jet\n");
    }

    //deposition in the coastline cell
    //for each grainsize-class add the sedflux to space
    int n = 0;
    double sedflux_sum = 0;
    for (n = 0; n < number_of_gscl; n++)
    {
      space[sim_time][i->row][i->col][n + 1] +=
        (sedflux_depo[j][n] * dt * dx);
      sedflux_sum += sedflux_depo[j][n];
    }
    space[sim_time][i->row][i->col][0] += (sedflux_sum * dt * dx);

  }

  for (; i != current_flowpath->end (); i++, j++)       // offshore in deltaic plume
  {

    printf ("=");

    double L_plume = (0.5 * b) + (tan_alpha * dx * (j - flx));
    int lat_width = (int) (L_plume / dy);

    //righthand side of the plume


    int o = 1;
    int graincell_cl = 0;
    for (int l = fly; l < (fly + (lat_width - 1)); l++, o++)
    {
      //lat_width is scaled to 0-2 domain appropriate for erf()
      l_scale[o] = (double) (l - fly) / (lat_width / 2.0);

      // the factor 0.5 needs to be used to divide half of the sediment only          
      double surf_fac_right = (erf (l_scale[o]) - erf (l_scale[o - 1])) * 0.5;
      //cout<<l<<" "<<"surf_fac_right = "<<surf_fac_right<<endl;

      //for each grainsize-class add the sedflux to space
      int n = 0;
      double sedflux_sum_r = 0;
      for (n = 0; n < number_of_gscl; n++)
      {
        sedflux_r[n] = 0.5 * sedflux_depo[j][n];
//                                      space[sim_time][i->row][l][n+1] +=(sedflux_depo[j][n]*surf_fac_right*dt*dx);
        sedflux_sum_r += sedflux_depo[j][n] * surf_fac_right * dt * dx;
      }

      space[sim_time][i->row][l][0] += sedflux_sum_r;

      // to divide the sediment from coarse (graincell_cl =0) to fine (graincell_cl = num_of_gscl)
      // assumption made: grainsizes-classes stack and don't mix
      double height = 0.0;
      while (height < sedflux_sum_r)
      {
        if (height + sedflux_r[graincell_cl] < sedflux_sum_r)

        {
          space[sim_time][i->row][l][graincell_cl + 1]
            += sedflux_r[graincell_cl];
          height += sedflux_r[graincell_cl++];
        }

        else
        {
          space[sim_time][i->row][l][graincell_cl + 1]
            += sedflux_sum_r - height;
          sedflux_r[graincell_cl] -= sedflux_sum_r - height;
          height = sedflux_sum_r;
        }
      }                         // end while
    }                           // end righthand

    //lefthandside of the plume
    int r = 1;
    graincell_cl = 0;
    for (int q = fly + 1; q > ((fly) - (lat_width - 2)); q--, r++)
    {
      double surf_fac_left = (erf (l_scale[r]) - erf (l_scale[r - 1])) * 0.5;
      //      cout<<q<<" "<<"surf_fac_left = "<<surf_fac_left<<endl;
      //for each grainsize-class add the sedflux to space

      int n = 0;
      double sedflux_sum_l = 0;

      for (n = 0; n < number_of_gscl; n++)
      {
        sedflux_l[n] = 0.5 * sedflux_depo[j][n];
        //      space[sim_time][i->row][q][n+1] += sedflux_depo[j][n]*surf_fac_left*dt*dx;
        sedflux_sum_l += sedflux_depo[j][n] * surf_fac_left * dt * dx;
      }

      //cout<<" j = "<<j<<" q ="<<q<<" sedfluxsum= "<<sedflux_sum<<endl;
      space[sim_time][i->row][q][0] += sedflux_sum_l;

      // to divide the sediment from coarse (graincell_cl =0) to fine (graincell_cl = num_of_gscl)
      // assumption made: grainsizes-classes stack and don't mix
      double height = 0.0;
      while (height < sedflux_sum_l)
      {
        if (height + sedflux_l[graincell_cl] < sedflux_sum_l)
        {
          space[sim_time][i->row][q][graincell_cl + 1]
            += sedflux_l[graincell_cl];
          height += sedflux_l[graincell_cl++];
        }
        else
        {
          space[sim_time][i->row][q][graincell_cl + 1]
            += sedflux_sum_l - height;
          sedflux_l[graincell_cl] -= sedflux_sum_l - height;
          height = sedflux_sum_l;
        }

      }                         //end while

    }                           //end righthand

  }                             // end of offshore flowpath loop

}                               //end of del_lat
