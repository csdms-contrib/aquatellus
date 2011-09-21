/*
 *  aquatellus_init_func.h
 *  
 *
 *  Created by Irina Overeem on 9/6/11.
 *  Copyright 2011 Univ of CO. All rights reserved.
 
 */

 #ifndef AQUATELLUS_INIT_FUNC_H_
 #define AQUATELLUS_INIT_FUNC_H_
 

 typedef double MATELEM;
 typedef MATELEM *vector;
 typedef vector *topo;
 typedef topo *grid;
 
 extern grid *space;

 void aquatellusreadinput();

 int define_space (int number_of_time_steps, int number_of_rows, int number_of_cols, int number_of_gscl);
 void initialize_space (double initial_height, double initial_gradient, double offshore_gradient);

 void aquatellus_init_func();

 
 #endif


