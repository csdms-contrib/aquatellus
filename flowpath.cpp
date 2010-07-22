#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "main.h"
#include "flowpath.h"
#include "Space.h"

using namespace std;

/****************************************************************
 SWITCH	: this function determines whether a switch in the river or delta
 flowpath occurs at that specific timestep. Sofar, this is based on tstep
 and a random factor;(lit-reference: Bahr et al, 2000) but it might be done 
 more process based in future. Whenever a switch occurs a new flowpath 
 is calculated, with the FLOWPATH module. 
 ****************************************************************/
double switch_flowpath()
{
	double noise2;
	noise2 = 1.0*rand()/RAND_MAX;//32767 is internaly defined for RAND_MAX
	int dt = 1;
	//cout<<noise2<<endl;
	// the time-step dependency is logarithmic
	double switch_threshold = log10(dt)/2.0 + noise2;
	cout << "switch_threshold" <<" "<<switch_threshold<< endl;
	return switch_threshold;
}

/****************************************************************************
 this determines which cell is the lowest cell in the adjacent cells to the 
 current cell. That cell is defined to be the cell to drain to
  **************************************************************************/

INTPAIR determine_lowest_cell (int t, INTPAIR actcel)
{
	INTPAIR lowercel;
	double max;
	double d = 0.0;
	int i,j, row_number, colum_number;
	
	row_number = actcel.row;
	colum_number = actcel.col;

	lowercel.row = -1;
	lowercel.col =-1;

	max = 0.0;

	// the entire matrix loop
	for (i=-1; i<=1; i++)
		{

			for (j=-1; j<=1;j++)	
				{
				if ((row_number+i>=0) && (row_number+i<number_of_rows) &&
					(colum_number+j>=0) && (colum_number+j<number_of_colums) &&
					!((i==0)&&(j==0))) 
				{
					d=get_topo_height (t,row_number, colum_number)
						- get_topo_height (t,row_number+i, colum_number+j);
					//cout << "Height difference is: " << d << endl;
					if((row_number!= row_number+i) && (colum_number!=colum_number+j))
					{
						d = d/sqrt(1.1);
					
						
						// this is used to give preference to diagonal drainage; 
						//if sqrt (1) is given; reaction is a more sinuous flowpath
						// if sqrt (2) is given; reaction is a more straight flowpath
					}
				}
				if(d>max)	
					{
						lowercel.row = row_number+i;
						lowercel.col = colum_number+j;
						max = d;
						//cout<<"max = "<<max<<endl;
					} 
				}//forloop
		}//forloop

	if ((lowercel.row ==-1) &&(lowercel.col ==-1))	
		{
		lowercel.row = actcel.row;
		lowercel.col = actcel.col;		
		}
	return lowercel;
}//end loop INTPAIR


/******************************************************************************
//the flowpath in the grid is determined based on the steepest descent algorithm
*********************************************************************************/

flowpath* determine_flowpath (int t)
{
	flowpath* new_flowpath = new flowpath ();
	INTPAIR current_cell, next_cell;



	current_cell.row = 0;
	// river starts in the top of the landscape in the middle cell
	current_cell.col = (int) (floor(number_of_colums/2));

	//river starts not every timestep in the same cell
	double randpos = 10.0*rand()/RAND_MAX; 
//	cout<<"randpos ="<<randpos<<endl;
	current_cell.col += (int)randpos;
//
	new_flowpath->push_back (current_cell);
	
//	cout << current_cell.row << " " << current_cell.col << endl;
	
	bool found = false;
	do	{
		next_cell = determine_lowest_cell(t, current_cell);

		// boundary condition for 3 landscape-edges
		// HIER NOG WEER NAAR KIJKEN (5 sept 2001)
		if ((next_cell.row == current_cell.row) &&
			(next_cell.col == current_cell.col) ||
				(next_cell.row == number_of_rows - 1) ||
			(next_cell.col == number_of_colums - 1) ||
			(next_cell.col == 0)) 
			{
				found = 1;
		} else {
			new_flowpath->push_back (next_cell);
		}

//		cout << next_cell.row << " "<< next_cell.col << endl;
		current_cell = next_cell;
	} while (!found);
	return new_flowpath;
}



flowpath* get_flowpath (int t)
{
	static flowpath* current_flowpath;

	if ((switch_flowpath () > 0.8)||(t==0))
		
	//if(t==0) 
	{
		delete current_flowpath;
		current_flowpath = determine_flowpath (t);
		cout <<"a new flowpath is determined" << endl;
	}
	else
	{
		cout << "no changes in the flowpath" << endl;
	}
	return current_flowpath;
}


