#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <math.h>
#include <time.h>
#include "aquatellus.h"
#include "main.h"

#define PI 3.14159265359 

//  RANDOM_CLIMATE FACTOR ////////////////////////////////////////////////////
// this factor induces changes over time in the initial input of discharge and
// initial sediment load --> reflecting a component of climatic changes in time
// It is set to random to facilitate research in which the climatic component
// is not that important 
///////////////////////////////////////////////////////////////////////////////////

double get_random_climate_factor()
{
	double upper_climate_range =1.6;
	climate_factor= upper_climate_range * rand()/RAND_MAX;			//random has range of 0-1.6
	//cout<<" climate-factor "<<climate_factor<<endl;
	return climate_factor;
}



double get_sealevel(int sim_time)
{
	///////////////////////////Sea level controls//////////////////////////////////
	double sea_lvl_ref=88;			// average level around which fluctuations occur
									// it is also the reference sealevel at t=0
	double sea_lvl_amp=4;			// amplitude in m
	double sea_lvl_prd = 800;		// period of the sine function in years

	sealevel = sea_lvl_ref+sea_lvl_amp*sin(2*PI*(sim_time/sea_lvl_prd)+PI);
	cout<< " sealevel = "<<sealevel<<endl;

	return sealevel;
}