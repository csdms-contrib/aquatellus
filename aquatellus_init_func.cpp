#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include "aquatellus.h"
#include "aquatellus_init_func.h"
#include "Space.h"
#include "main.h"

using namespace std;

void
aquatellus_init_func ()
{

  //READ INPUT

  aquatellusreadinput ();

  // INITIALIZE THE LANDSCAPE AND SEDIMENT GRID

  // step 1 initialize the memory for the grid, and the grainsize classes that fill it over time
  if (define_space
      (number_of_time_steps, number_of_rows, number_of_colums,
       number_of_gscl) == 0)
  {
    cout << "space has been declared!" << endl;
  }

  else
  {
    cout << "Careful! not enough memory for space_array" << endl;
  }

  // step 2 initialize the actual values that fill the SEDIMENT GRID and define the topography/bathymetry
  initialize_space (initial_height, initial_gradient, offshore_gradient);

}
