#ifndef OUTPUT_H_
#define OUTPUT_H_

// min -max comparisons needed for write grainsize plot
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define EPSloc 0.00001


void open_data_files(FILE **ftopo, FILE **fpath, int time = 0);
void close_data_files(FILE **ftopo, FILE **fpath);
void write_topographical_map_header(FILE **ftopo);
void write_topographical_map(FILE **ftopo, long sim_time);
double ** calculate_grain_distr(double * grain_size, int Xpos);
double * calculate_grain_distr_well(double *grain_size, int wellpos_row, int wellpos_col);
void write_median_psfile(double **median, int Xpos);
void write_simulatedwell_psfile(double *median_point, int wellpos_row,int wellpos_col);

#endif
