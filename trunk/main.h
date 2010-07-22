#ifndef _MAIN_H_
#define _MAIN_H_

//Here the globals of the program are declared

// global that initialize the simulation dimensions
extern int number_of_time_steps;
extern int number_of_rows;
extern int number_of_colums;
extern int number_of_gscl;
extern int dx;
extern int dy;
extern int dt;

extern double sed_cont_pct[];

extern long sim_time;
extern long end_of_times;

extern double sealevel;
extern double climate_factor;

// globals for erosion and sedimentation calculations
extern int coastline_cell;

#endif
