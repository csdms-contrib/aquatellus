/*
 *  run_aquatellus.cpp
 *
 *
 *  Created by Irina Overeem on 9/9/11.
 *  Copyright 2011 Univ of CO. All rights reserved.
 *
 */

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
#include "aquatellus_finalize.h"
#include "run_aquatellus.h"

using namespace std;

FILE *ftopo;

FILE *fpath;

void
run_aquatellus ()
{
  printf ("simtime=%d\n", sim_time);
  // FILE I/O
  open_data_files (&ftopo, &fpath, sim_time);

  //  CONTROLLING FACTORS

  // controlling time-dependent factors are determined
  double sealevel = get_sealevel (sim_time);

  get_random_climate_factor ();

  // SWITCH + FLOWPATH is determined
  flowpath *current_flowpath = get_flowpath (sim_time);

  flowpath::iterator i;

  // current flowpath is written to file and array
  int k = 0;

  for (i = current_flowpath->begin (); i != current_flowpath->end (); i++, k++)
  {
    // write flowpath to array
    //                      flowpath_array[k][0]=i->row;
    //                      flowpath_array[k][1]=i->col;

    //cout<<flowpath_array[k][0]<<" "<<flowpath_array[k][1]<<endl;
    //write flowpath to file
    fprintf (fpath, "%d %d\n", i->row, i->col);
    //      cout<<sim_time<<" "<<i->row<<" "<<i->col<<endl;

  }

  // 2D DYNAMICS ALONG THE LONGITUDINAL PROFILE
  int j;

  // coastline cell is determined
  flowpath::iterator coastline_cell_index =
    determine_coastline_cell (current_flowpath, sealevel, sim_time);
  cout << "Current coastline cell index has row: " << coastline_cell_index->row
    << " and column " << coastline_cell_index->col << endl;

  // discharge
  double *discharge = determine_discharge_along_flowpath (current_flowpath,
                                                          coastline_cell_index,
                                                          sealevel,
                                                          sim_time);

  /*      cout << " Found the following discharge array along flowpath: " << endl;
     cout << discharge[0];
     for (j = 1; j < current_flowpath->size (); j++) { //  Only for debugging purposes
     cout << ", " << discharge[j];
     }
     cout << endl; */

  // erosion
  double *dh_erosion = determine_erosion_rate_along_flowpath (current_flowpath,
                                                              coastline_cell_index,
                                                              discharge,
                                                              sim_time);

  //cout << " Found the following erosion heights along flowpath: " << endl;
  //cout << dh_erosion[0];
  //for (j = 1; j < current_flowpath->size (); j++) { //  Only for debugging purposes
  //      cout << ", " << dh_erosion[j];
  //}
  //      cout << endl;

  FILE *afile = fopen ("eroflux.txt", "w");

  double **sedflux_ero =
    dig_erosion_from_space (current_flowpath, dh_erosion, sim_time);
  fprintf (afile,
           "Dug the following amounts of sediment from space along flowpath: \n");
  for (j = 0; j < current_flowpath->size (); j++)
  {     //  Only for debugging purposes
    for (int k = 0; k < number_of_gscl; k++)
    {
      fprintf (afile, ",%f", sedflux_ero[j][k]);
    }
    fprintf (afile, "\n");
  }
  fclose (afile);

  // sedimentation flux
  FILE *bfile = fopen ("finalflux.txt", "w");

  double **sedflux_depo = sedimentation_rate_along_flowpath (current_flowpath,
                                                             coastline_cell_index,
                                                             sedflux_ero);
  fprintf (bfile, "Final sediment fluxes along flowpath: \n");
  for (j = 0; j < current_flowpath->size (); j++)
  {     //  Only for debugging purposes
    for (int k = 0; k < number_of_gscl; k++)
    {
      fprintf (bfile, ",%f", sedflux_depo[j][k]);
    }
    fprintf (bfile, "\n");
  }
  fclose (bfile);

  // lateral distribution of sedflux in fluvial domain
  //      double *firstder = calculate_firstder(current_flowpath);
  //      double *y2 = spline(current_flowpath, firstder);
  //      double *curvature =calculate_curvature(current_flowpath, firstder, y2);

  //      double *asymmetry = calculate_asymmetry(current_flowpath, curvature);
  ///     for (j = 0; j < current_flowpath->size (); j++) { //  Only for debugging purposes
  ///     cout << ","<<asymmetry [j];
  //      }
  //      cout << endl;

  fluv_lat (current_flowpath, coastline_cell_index, sim_time, discharge,
            sedflux_depo);

  //lateral distribution of sediment in the deltaic domain
  del_lat (current_flowpath, coastline_cell_index, discharge, sim_time,
           sedflux_depo);

  write_topographical_map (&ftopo, sim_time);

  free (discharge);
  free (dh_erosion);
  free (sedflux_ero);
  free (sedflux_depo);
  //      free(firstder);
  //      free(y2);
  //      free(curvature);
  //      free(asymmetry);

  close_data_files (&ftopo, &fpath);

}
