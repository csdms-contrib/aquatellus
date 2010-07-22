#ifndef SPACE_H_
#define SPACE_H_

///////////////////////////////////////////////////////////////////////////
// Define Space
// 
// 	this function generates a landscape as a matrix with a specified
// number of rows and colums. The general trend is a decreasing gradient
// (dy) towards the sea and shelf edge. To incorporate the local variations
// in topographical height a small random height is added to the generated
// trend topography to every grid cell
///////////////////////////////////////////////////////////////////////////

typedef double MATELEM;
typedef MATELEM *vector;
typedef vector *topo;
typedef topo *grid;

extern grid *space; 
int define_space(int number_of_time_steps,
					 int number_of_rows,
					 int number_of_cols,
					 int number_of_gscl);
void initialize_space();
double get_space_height(int row, int column);
double get_space_height (int t, int row, int column);
double get_topo_height(int t, int row, int col);
double grain_percent(int t, int row, int column, int grainclass);

#endif
