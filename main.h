#ifndef MAIN_H_
#define MAIN_H_

//Here the globals of the program are declared

// global that initialize the simulation dimensions
extern int number_of_time_steps;

extern int dt;

extern long sim_time;

extern long end_of_times;

// GRID GLOBALS
extern int number_of_rows;

extern int number_of_colums;

extern int number_of_gscl;

extern int dx;

extern int dy;

extern double initial_height;

extern double initial_gradient;

extern double offshore_gradient;

//SEDIMENT GLOBALS
extern double grain_size[];

extern double sed_cont_pct[];

extern double k_er_fluv;

extern double k_er_marine;

extern double traveldist_fluvial[];

extern double traveldist_marine[];

//CONTROL FACTORS GLOBALS
extern double discharge_def;

extern double sediment_load_def;

extern double sealevel;

extern double sea_lvl_ref;

extern double sea_lvl_amp;

extern double sea_lvl_prd;

extern double climate_factor;

// globals for erosion and sedimentation calculations
extern int coastline_cell;

extern int Xpos;

extern int wellpos_row;

extern int wellpos_col;

#endif
