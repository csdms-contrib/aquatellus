#ifndef DEL_LAT_H_
#define DEL_LAT_H_

void del_lat (flowpath * current_flowpath,
              flowpath::iterator coastline_cell_index, double *discharge,
              long sim_time, double **sedflux);

#endif
