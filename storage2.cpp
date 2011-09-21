/********************************************************************

  STORAGE this function calculates the storage capacity of the
  cross sectional profile based on its topography and the cross
  sectional discharge .

  its only a test version for incorporation in FLUV_LAT
*********************************************************************/
/*
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
# include <time.h>
#include <iterator>
#include "main.h"
#include "aquatellus.h"
#include "space.h"
#include "flowpath.h"

typedef struct {
		int left_count;
		int right_count;
	} COORPAIR;

COORPAIR * storage(flowpath * current_flowpath, double * discharge)

{
	int k;
	int fly;
	int right_count,left_count, lat_width_right,lat_width_left;

	
	double A_cum; // the cumulative cross sectional discharge that can be stored
	double waterlevel_height, old_waterlevel_height;

	

	COORPAIR curve_width;

	//double new_waterlh =0; // the waterlevel height after fill
	//double right_bankh=0; // for right handside river bank height
	//double left_bankh =0; // for left handside river bank height

	int j=0;
	for (flowpath::iterator i = current_flowpath->begin(); i != current_flowpath->end (); i++, j++)
	 {

		//initialize flowpath position
		
		left_count = right_count = i->col;
		waterlevel_height= get_space_height(sim_time,i->row, i->col); // bottom of the river is the initial water level height
		A_cum =0;

	
	// loop runs as long as the stored discharge is less than the total cross sectional Q
	while (A_cum< discharge[j])
	{
			// determine the riverbed
		    while (get_space_height(sim_time, i->row, right_count) <= waterlevel_height)
		    right_count++;
					
			while (get_space_height(sim_time, i->row, left_count) <= waterlevel_height)
			left_count--;
			
			// determine whether rightside is the lowest riverbank					
			if (get_space_height(sim_time, i->row,right_count) <= get_space_height(sim_time, i->row,left_count))
				{
					waterlevel_height = get_space_height (sim_time, i->row, right_count);

					// from right riverbank onwards the cross-sectional storage is
					// calculated untill the leftbank is reached
					for (int k = right_count-1; k>left_count; k --)
					{A_cum += ((waterlevel_height- get_space_height(sim_time, i->row,k))* dy);}

				} //endif

					// otherwise leftside is the lowest riverbank
			        else	{
							waterlevel_height = get_space_height(sim_time, i->row, left_count);
							// in that case from the left riverbank onwards the
							// cross sectional Q is determined
							for (int k = left_count+1; k< right_count; k ++)
							{A_cum += ((waterlevel_height- get_space_height(sim_time, i->row,k))* dy);}
					
					}
						
 				
					if (A_cum > discharge[j])  break;
					
					// each step all area under the waterlevel unto the topographical
					// height is calculated, thus the A_cum has to be set back to 0

					A_cum =0;

				

	}// end of while loop

	
        lat_width_left = fly - left_count;
		lat_width_right = right_count -fly;
	//	cout <<lat_width_left<< " "<<  lat_width_right<< endl;

		
		//printf(" %d %d \n", curve_width.left_count, curve_width.right_count);
        return curve_width;

	}

}*/
