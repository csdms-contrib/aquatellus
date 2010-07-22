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
#include "Space.h"
#include "flowpath.h"
#include "inputcontrols.h"
#include "dynamics.h"
#include "fluv_lat.h"
#include "del_lat.h"
#include "output.h"
#include "ei.h"

using namespace std;

// time & grid dimensions
// it is recommended to keep them rather low to avoid long simulation runs

long end_of_times = 200;
int number_of_rows = 120;
int number_of_colums = 120;
int number_of_gscl = 4;
int dx =1;
int dy =1;
int dt = 1;
int number_of_time_steps = end_of_times/dt;

double grain_size[]= {0.400,0.200,0.012,0.004}; // grainsizes in mm, coarse to fine
double sed_cont_pct[]={0.25,0.25,0.25,0.25};	//percentage of each grainsize
int Xpos =85;				 //position at which Xsection for output is made
int num_outp = 20;			// number of times that topography file gets written away
// X-Y position of simulated wells 
int wellpos_row =45;	 
int wellpos_col =85;

double *discharge;
double *dh_erosion;
int **flowpath_array;
double *firstder;
COORPAIR curve_width;

long sim_time;
double sealevel;
double climate_factor;

FILE *ftopo;
FILE *fpath;


///////////////////// THE ACTUAL SIMULATION//////////////////////////////////
int main ()

{	
	

	// LANDSCAPE INITIALIZATION
	if (define_space (number_of_time_steps, number_of_rows, number_of_colums,number_of_gscl ) == 0)
		
	{
		cout << "space has been declared!" << endl;
	}
	else
	{
		cout << "Careful! not enough memory for space_array" << endl;
	}

	initialize_space();

	

long sim_time = 0;

//time loop
while (sim_time<end_of_times) 
	{
	
	
	// FILE I/O
	open_data_files(&ftopo, &fpath,sim_time);
	write_topographical_map_header(&ftopo);
	//	INPUT DATA CONTROLLING FACTORS
	// controlling time-dependent factors are determined
	double sealevel = get_sealevel(sim_time);
	get_random_climate_factor();
	//cout << climate_factor<< endl;	

	// SWITCH + FLOWPATH is determined
	flowpath* current_flowpath = get_flowpath (sim_time);
	flowpath::iterator i;

	//dynamically declare flowpath-array
	// h gives the row and col dimension (2 array elements)
	//int h=2;
	//flowpath_array = (int**)calloc(current_flowpath->size(), sizeof(int));
//	if(flowpath_array==NULL)
//		return -1;
//
//		for (int m=0;m<current_flowpath->size();m++)//
//
//		{
//			flowpath_array[m]=(int*)calloc(h,sizeof(int));
//			if(flowpath_array ==NULL)
//				return -1;
//		}

	// current flowpath is written to file and array
	int k =0;
	for (i = current_flowpath->begin (); i != current_flowpath->end (); i++, k++)
		{
			// write flowpath to array
//			flowpath_array[k][0]=i->row;
//			flowpath_array[k][1]=i->col;
			
			//cout<<flowpath_array[k][0]<<" "<<flowpath_array[k][1]<<endl;
			//write flowpath to file
			fprintf(fpath, "%d %d\n", i->row, i->col);
		//	cout<<sim_time<<" "<<i->row<<" "<<i->col<<endl;
			
		}
		


	// 2D DYNAMICS ALONG THE LONGITUDINAL PROFILE 
	int j;
	// coastline cell is determined
	flowpath::iterator coastline_cell_index = determine_coastline_cell (current_flowpath, sealevel, sim_time);
	cout << "Current coastline cell index has row: " 
		<< coastline_cell_index->row
		<< " and column " 
		<< coastline_cell_index->col 
		<< endl;

	// discharge
	double* discharge = determine_discharge_along_flowpath (current_flowpath, coastline_cell_index, sealevel, sim_time);
/*	cout << " Found the following discharge array along flowpath: " << endl;
	cout << discharge[0];
	for (j = 1; j < current_flowpath->size (); j++) { //  Only for debugging purposes
		cout << ", " << discharge[j];
	}
	cout << endl;*/
	
	// erosion
	double* dh_erosion = determine_erosion_rate_along_flowpath(current_flowpath, coastline_cell_index, discharge, sim_time);
	//cout << " Found the following erosion heights along flowpath: " << endl;
	//cout << dh_erosion[0];
	//for (j = 1; j < current_flowpath->size (); j++) { //  Only for debugging purposes
	//	cout << ", " << dh_erosion[j];
	//}
//	cout << endl;

	FILE* afile = fopen ("eroflux.txt", "w");
	double** sedflux_ero = dig_erosion_from_space(current_flowpath, dh_erosion, sim_time);
	fprintf (afile, "Dug the following amounts of sediment from space along flowpath: \n");
	for (j = 0; j < current_flowpath->size (); j++) { //  Only for debugging purposes
			for (int k = 0; k <number_of_gscl; k++) {
			fprintf (afile, ",%f", sedflux_ero[j][k]);
		}
		fprintf (afile, "\n");
	}
	fclose (afile);

	// sedimentation flux
	FILE* bfile = fopen ("finalflux.txt", "w");
	double** sedflux_depo = sedimentation_rate_along_flowpath(current_flowpath, coastline_cell_index, sedflux_ero);
	fprintf (bfile, "Final sediment fluxes along flowpath: \n");
	for (j = 0; j < current_flowpath->size (); j++) { //  Only for debugging purposes
			for (int k = 0; k <number_of_gscl; k++) {
			fprintf (bfile, ",%f", sedflux_depo[j][k]);
		}
		fprintf (bfile, "\n");
	}
	fclose (bfile);
	
	// lateral distribution of sedflux in fluvial domain
//	double *firstder = calculate_firstder(current_flowpath);
//	double *y2 = spline(current_flowpath, firstder);
//	double *curvature =calculate_curvature(current_flowpath, firstder, y2);

	
//	double *asymmetry = calculate_asymmetry(current_flowpath, curvature);
///	for (j = 0; j < current_flowpath->size (); j++) { //  Only for debugging purposes
///	cout << ","<<asymmetry [j];
//	}
//	cout << endl;

	fluv_lat(current_flowpath,coastline_cell_index, sim_time, discharge, sedflux_depo);

	//lateral distribution of sediment in the deltaic domain
	del_lat(current_flowpath,coastline_cell_index, discharge,sim_time, sedflux_depo);

	
	write_topographical_map(&ftopo, sim_time);


	free (discharge);
	free (dh_erosion);
	free (sedflux_ero);
	free (sedflux_depo);
//	free(firstder);
//	free(y2);
//	free(curvature);
//	free(asymmetry);

	
	if (sim_time>=end_of_times)
	{
		printf("end of programme\n");
	}

	sim_time+=dt;

	close_data_files(&ftopo, &fpath);
	} //end of simulation time loop 
	
	//for Xsection output
	// X-section output file as post-processing (this can be made much more sophisticated)
	double ** median = calculate_grain_distr(grain_size, Xpos);
	// for well output
//	double * median_point =calculate_grain_distr_well(grain_size, wellpos_row,wellpos_col);

	//	//for Xsection output
	write_median_psfile(median, Xpos);
		//for well output
//	write_simulatedwell_psfile(median_point, wellpos_row,wellpos_col);

	
	return 0;

}

