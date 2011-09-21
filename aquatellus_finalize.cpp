/*
 *  aquatellus_finalize.cpp
 *
 *
 *  Created by Irina Overeem on 9/8/11.
 *  Copyright 2011 Univ of CO. All rights reserved.
 *
 */

#include "aquatellus_finalize.h"
#include "aquatellus.h"
#include "output.h"
#include "main.h"

void
aquatellus_finalize ()
{
  // POST-PROCESSING OF THE DATA CUBE

  //for Xsection output
  // X-section output file as post-processing
  double **median = calculate_grain_distr (grain_size, Xpos);

  // for well output
  double *median_point =
    calculate_grain_distr_well (grain_size, wellpos_row, wellpos_col);

  // WRITE THE OUTPUT FILES

  //for Xsection output
  write_median_psfile (median, Xpos);

  //for well output
  write_simulatedwell_psfile (median_point, wellpos_row, wellpos_col);

}
