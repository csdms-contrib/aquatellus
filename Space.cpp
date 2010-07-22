#include <malloc.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include "main.h"
#include "Space.h"
#include "aquatellus.h"

grid *space;


//////////////////////////////////////////////////////////////////////////////
// here the 4-D matrix (time and 3 spatial dimensions) is declared
//////////////////////////////////////////////////////////////////////////////

int define_space(int number_of_time_steps,
					 int number_of_rows,
					 int number_of_colums,
					 int number_of_gscl)
					 
{
	space = (grid*)calloc (number_of_time_steps, sizeof (grid));
	if (space == NULL)
		return -1;

	for (int t = 0; t < number_of_time_steps; t++)
	{
		space[t] = (grid)calloc (number_of_rows, sizeof (topo));
		if (space[t] == NULL)
			return -1;

		for (int i = 0; i < number_of_rows; i++)
		{
			space[t][i] = (topo) calloc (number_of_colums, sizeof (vector));
			if (space[t][i] == NULL)
				return -1;

			for (int j = 0; j < number_of_colums; j++)
			{
				space[t][i][j] = (vector) calloc (number_of_gscl+1, sizeof (MATELEM));
				if (space[t][i][j] == NULL)
					return -1;
			}
		}
	}
	return 0;
}

///////////////////////////////////////////////////////////////////////////
/* initialize space

This function generates a topographical space as a matrix with a specified
number of rows and colums. The general trend is a decreasing gradient
(dy) towards the sea and shelf edge. 
To incorporate the local variations in topographical height a small random
height is added to the generated trend topography to every grid cell
Subsequently,the initial height at each x,y position is filled with sediments 
of different grainsize classes. The dh distribution depends on user-defined
percentages (sed_cont_perc)
*/
///////////////////////////////////////////////////////////////////////////

void initialize_space()
{ 
	int i,j,k;
	double noise;
	double initial_gradient = 0.2;	// initial gradient in longitudinal direction
	double offshore_gradient= 0.2;
	
	// this adds a random small height to each x,y topo point of the matrix
	/* Seed the random-number generator with current time so that
	 the numbers will be different every time we run. */   
//	srand ((unsigned)time (NULL));

	for (j = 0; j < number_of_colums; j++)
	{
		space[0][0][j][0]=100.0; // topographical height of top of the profile
		noise = rand ()*(initial_gradient/2.0);
		noise = noise/100000.0;
		space[0][0][j][0] += noise;
	}

		for (i=1; i<number_of_rows-30; i++)
	{
		for(j=0;j<number_of_colums;j++)
		{
			space[0][i][j][0] = space[0][i-1][j][0] - initial_gradient;

			//cout<< space[0][i][j][0]<< endl;	
				
			noise = rand ()*(initial_gradient/2.0);
			noise = noise/100000.0;

		
			space[0][i][j][0]=space[0][i][j][0]+noise;
			// check 1 total heights's declared?
			//cout<< space[0][i][j][0]<< endl;

			// fill with dh's for each grainsize class
			for (k=1; k<=number_of_gscl; k++)
			{
				space[0][i][j][k] = space[0][i][j][0]*sed_cont_pct[k-1];
				// check 2 partial heights declared?
				//cout<< space[0][i][j][k]<< endl;
			}
		}
	}

		for (i=number_of_rows -30; i<number_of_rows; i++)
	{
		for(j=0;j<number_of_colums;j++)
		{
			space[0][i][j][0] = space[0][i-1][j][0] - offshore_gradient;

			//cout<< space[0][i][j][0]<< endl;	
				
			noise = rand ()*(initial_gradient/2.0);
			noise = noise/100000.0;
			space[0][i][j][0]=space[0][i][j][0]+noise;
			// check 1 total heights's declared?
			//cout<< space[0][i][j][0]<< endl;

			// fill with dh's for each grainsize class
			for (k=1; k<=number_of_gscl; k++)
			{
				space[0][i][j][k] = space[0][i][j][0]*sed_cont_pct[k-1];
				// check 2 partial heights declared?
				//cout<< space[0][i][j][k]<< endl;
			}
		}
	}
}

double get_space_height(int row, int column)
{
	return space[0][row][column][0];
}

// get space height of previous time step.
double get_space_height (int t, int row, int column)
{
	return space[t][row][column][0];
}

// the actual topographic height is necessary to calculate the slopes
// in the erosion flux calculation, to compare with waterdepth etc
double get_topo_height(int t, int row, int col)
{
	int k=0;
	double topoheight=0;
	for (k=0; k<=t; k++)
	{
		topoheight += space[k][row][col][0];
	}
	
	return topoheight;
}

double grain_percent(int t, int row, int column, int grainclass)
{
	return space[t][row][column][grainclass]/space[t][row][column][0];
}
