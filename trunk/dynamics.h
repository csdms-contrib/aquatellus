#ifndef DYNAMICS_H_
#define DYNAMICS_H_

flowpath::iterator determine_coastline_cell(flowpath* current_flowpath, double sealevel, long sim_time);
double* determine_discharge_along_flowpath(flowpath* current_flowpath, flowpath::iterator coastline_cell_index, double sealevel, long sim_time);
double* determine_erosion_rate_along_flowpath(flowpath* current_flowpath, flowpath::iterator coastline_cell_index, double* discharge, long sim_time);
double** dig_erosion_from_space(flowpath* current_flowpath, double* dh_erosion, long sim_time);
double ** sedimentation_rate_along_flowpath(flowpath* current_flowpath, flowpath::iterator coastline_cell_index, double** sedflux_ero);

#endif

