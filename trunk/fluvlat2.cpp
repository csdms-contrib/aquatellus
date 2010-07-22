///////////////////////////////////////////////////////////////////////////////////////////
//FLUV-LAT (FLUVial sedimentation distribution in LATeral direction)
//
//-the first and second derivatives are calculated and the spline
// for an file with x,y data. This file is the flowpath of the river in 
// the landscape matrix.
//
//- the curvature of the splined flowpath is calculated to asses
//  the sinuosity of the channel belt. 
// 
//-based on the curvature asymmetry in sediment distribution 
// over cross profiles is calculated.
//
//-based on storage capacity of the local x-profile topography and
// the x-sectional discharge; the lateral width of the sedimentation 
// curve is determined.
// 
//- the EI, exponential integral, determines the relative areas under
// the sedimentation curve.
//
//- nett lateral sedimentation is calculated.
//
///////////////////////////////////////////////////////////////////////////////////////////
#include "fluv_lat.h"
#include "Nr.h"
#include "Nrutil.h"


////////////////////////////////////////////////////////////////////////
// this function calculates the first derivatives of the x,y data that form
// the flowpath.Code is based on Stirlings formula (averages derivative calculation
// in forward and backward direction). Lit ref: "Numerical Recipes in C"
// (Press et al, 1996) and "C programming; the essentials for engineers and 
// scientists" (Brooks, 1998).
// The first derivatives of flowpath are needed for the calculation of the curvature;
///////////////////////////////////////////////////////////////////////

double* calculate_firstder(flowpath* current_flowpath)
{
	double* firstder;
	if ((firstder=(double*)calloc(current_flowpath->size (),sizeof(double)))==NULL)
		printf("error in initializing firstder array");

	int j=0;
	for (flowpath::iterator i = current_flowpath->begin(); i != current_flowpath->end (); i++, j++) {	
		if (j == 0) {
			// boundary condition for begin of input range, when the step backward
			//values of ya are not known 
			firstder[j]=0.0;
		} else if (i != current_flowpath->end()) {
			// boundary condition for end of input range, when the step forward
			//values of ya are not known anymore
			flowpath::iterator prev = i;
			prev--;
			flowpath::iterator next = i;
			next++;
			if (next != current_flowpath->end ()) 
			{
				firstder [j] = ((i->col - prev->col) / (i->row - prev->row) + (next->col - i->col) / (next->row - i->row)) / 2;
			} 
			else 
			{
				if ((i->row == prev->row)||(next->row == i->row))
					firstder[j] = firstder[j - 1]; // Not strictly the first derivative
				else 
					firstder[j] = ((i->col - prev->col) / (i->row - prev->row)) / 2;
			}
		}
			
	}

	return firstder;
}




//////////////////////////////////////////////////////////////////////////////////////////
// SPLINE ROUTINE 1
// given arrays x[1...n] and y[1...n] containing a tabulated function, i.e. yi = f(xi),
//with x1<x2<x3...<xn and given values yp1 and ypn for the first derivative of the 
//interpolating function at point 1 and n; this routine returns an array y2 [1..n] that
//contains the second derivative of the interpolating function at the tabulated points xi.
//If yp1 and/or ypn are equal 1*10 ^ 30 or larger, the routine is signaled to set the corresponding
//boundary condition for a natural spline, with zero derivative on that boundary
///////////////////////////////////////////////////////////////////////////////////////////

double * spline(flowpath* current_flowpath, double* firstder)
{
	int k, np;
	double p, qn, sig, un, *u;
	np = current_flowpath->size();
	u = new double [np];

	double* y2;
	if ((y2=(double*)calloc(current_flowpath->size (),sizeof(double)))==NULL)
		printf("error in initializing secder array");


	flowpath::iterator i = current_flowpath->begin ();
	int j = 0;
	if (firstder[0] > 0.99e30)
		{
			y2[0]=u[0]=0.0;
			// the lower boundary is set to natural
			// or else to specified first derivative
		}
	else
		{
			flowpath::iterator next = i;
			next++;

			y2[0] = -0.5;
			
			if (i->row == next->row)
			
				u[0] = (3.0 /(0.5)) *((next->col - i->col)/(0.5)-firstder[j]);
			else
				u[0] = (3.0 /(next->row - i->row)) *((next->col - i->col)/(next->row-i->row)-firstder[j]);
		}

	// this is the decompositional loop of the tridiagonal algorithm. y2 and u
	// are used for temporary storage of the decomposed factors
	for (i++,j++; i != current_flowpath->end (); i++, j++)
	{
		flowpath::iterator prev = i;
		prev--;
		flowpath::iterator next = i;
		next++;

		if (i != current_flowpath->end()) 
			{
		
				if ((i->col == prev->col) || (i->col ==next->col)||(next->row==prev->col))
				{
					
					y2[j] = 1e029;
				}

				else
				{
					sig = (i->row - prev->row) / (next->row - prev->col); 
					p = sig * y2[j-1] + 2.0;
					y2[j] = (sig - 1.0) /p;
					
					if ((i->row == next->row)||(i->row == prev->row))
						u[j] = ( next->col - i->col)/ (0.5) - (i->col - prev->col)/ (0.5);

					else
						u[j] = (( next->col - i->col)/ (next->row - i->row)) - ((i->col - prev->col)/ (i->row - prev->row));
				
					u[j] = (6.0*u[j]) / (next->row - prev->row) - sig*u[j-1]/p;
				}

			}


		else 
			{

			flowpath::iterator prev2 = prev;
			prev2--;

			if (firstder[np] >0.99e30)
				qn=un=0.0;
				else 
					{
						qn=0.5;
						u[np-1] = (3.0/(prev->row - prev2->row))*(firstder[np]-(prev->col -prev2->col)/(prev->row-prev2->row));
					}
				y2[np-1] = (un-qn*u[np-2]) / (qn*y2[np-2] + 1.0);
				for (k=np-2; k>=0; k--)
				y2[k] = y2[k]*y2[k+1]+u[k];

			}

	

	}

		delete [] u;
	return y2;
}



/////////////////////////////////////////////////////////////////////////////////////////
// This function calculates the curvature of the flowpath. The function reads the first
// (first-der) and second derivatives (sec-der)of the flowpath function as previously 
// calculated in SPLINE, which can then be used to calculate the curvature. The formula 
// is retrieved from:Tallarida, R.J. (1992) Pocketbook of integrals and mathematical formulas
// 2 nd edition.
/////////////////////////////////////////////////////////////////////////////////////////

double * calculate_curvature(flowpath * current_flowpath, double * firstder, double *y2)

{
	double* curvature;
	if ((curvature=(double*)calloc(current_flowpath->size ()+1,sizeof(double)))==NULL)
		printf("error in initializing curvature array");
	
	int j=0;
	for (flowpath::iterator i = current_flowpath->begin(); i != current_flowpath->end (); i++, j++)
	{
		// calculation of curvature
		firstder[j]= pow(firstder[j],2.0);
		firstder[j] = 1+firstder[j];
		curvature[j] = y2[j]/(pow(firstder[j],(3./2.)));
	}

  return curvature;
}


////////////////////////////////////////////////////////////////////////
// In this function the curvature is normalized and the asymmetry 
// coefficient determined
///////////////////////////////////////////////////////////////////////

double * calculate_asymmetry(flowpath * current_flowpath, double * curvature)

{

	double* asymmetry;
	if ((asymmetry=(double*)calloc(current_flowpath->size (),sizeof(double)))==NULL)
		printf("error in initializing asymmetry array");
	

	double min_curvature = 0;
	double max_curvature = 0;
	double range =0;	

	int j=0;

	for (flowpath::iterator i = current_flowpath->begin(); i != current_flowpath->end (); i++, j++)
		{
	
			if (curvature[j]> max_curvature)
				{max_curvature = curvature[j];}

			if (curvature[j]< min_curvature)
				{min_curvature = curvature[j];}
				//cout<<max_curvature<<" "<<min_curvature<<endl;
		}

		range = (max_curvature - min_curvature);

	j=0;
	for (flowpath::iterator i = current_flowpath->begin(); i != current_flowpath->end (); i++, j++) 
		{

			asymmetry[j] = ((curvature[j] - max_curvature)+(curvature[j]-min_curvature))/range;

		}

		return asymmetry;
}
 

/********************************************************************

  STORAGE this function calculates the storage capacity of the
  cross sectional profile based on its topography and the cross
  sectional discharge .
*********************************************************************/

COORPAIR storage(flowpath::iterator i, double * discharge, int j, long sim_time)

{
	
		int right_count,left_count, lat_width_right,lat_width_left;
		double A_cum; // the cumulative cross sectional discharge that can be stored
		double waterlevel_height;
			
		COORPAIR curve_width;		
		
			
		// initialize flowpath position
		left_count = right_count = i->col;
		waterlevel_height= get_topo_height(sim_time,i->row, i->col); // bottom of the river is the initial water level height
		A_cum =0;
					
		// loop runs as long as the stored discharge is less than the total cross sectional Q
		// and boundary condition for if the entire matrix width is flooded
		// especially in the beginning; height differences will be small
		// not much Q can be stored.

		while (A_cum<discharge[j])
			{
				// determine the riverbed
				while ((get_topo_height(sim_time, i->row, right_count) <= waterlevel_height)&& (right_count< number_of_colums-1))
					right_count++;
				
				while ((get_topo_height(sim_time, i->row, left_count) <= waterlevel_height) &&(left_count >0))
					left_count--;
				
				// determine whether rightside is the lowest riverbank					
				if (get_topo_height(sim_time, i->row,right_count) <= get_topo_height(sim_time, i->row,left_count)) 
				{
					waterlevel_height = get_topo_height (sim_time, i->row, right_count);
					
					// from right riverbank onwards the cross-sectional storage is
					// calculated untill the leftbank is reached
					for (int k = right_count-1; k>left_count; k --)
					{
						A_cum += ((waterlevel_height- get_topo_height(sim_time, i->row,k))* dy);
					} 
					
				} //endif
				
				// otherwise leftside is the lowest riverbank
				else	
				{
					waterlevel_height = get_topo_height(sim_time, i->row, left_count);
					// in that case from the left riverbank onwards the
					// cross sectional Q is determined
					for (int k = left_count+1; k< right_count; k ++)
					{
						A_cum += ((waterlevel_height- get_topo_height(sim_time, i->row,k))* dy);
					}
					
				}
				
				// lat width depends on the storage capacity
				if ((A_cum > discharge[j]) || (right_count==number_of_colums-1)||(left_count ==0))
				{ 
					break; 
				}
			
				// for each step along the flowpath 
				// the area under the waterlevel unto the topographical
				// height is calculated, thus the A_cum has to be set back to 0
				A_cum =0;
							
			}// end of while loop
		
			lat_width_left = i->col - left_count;
			lat_width_right = right_count -i->col;
			
		//	curve_width.left = lat_width_left;
		//	curve_width.right= lat_width_right;
			curve_width.left =left_count;
			curve_width.right= right_count;
		
	return curve_width;
}	
// Dit is de sept 2001 versie waarin 'storage gebruikt wordt om de depositie zones te bepalen; er lijken nog wel wat bugs in te zittten
// Daarom maar even voor de simpelere bovenstaande oplossing gekozen.
/*void fluv_lat(flowpath* current_flowpath,flowpath::iterator coastline_cell_index, long sim_time, COORPAIR curve_width, double *discharge, double **sedflux_depo)

{
	double *cellpartright;
	cellpartright = (double *) calloc(number_of_colums,sizeof(double));

	double *cellpartleft;
	cellpartleft = (double *) calloc(number_of_colums,sizeof(double));

		// to divide the depositional sediment flux over a left and righthand component
	double *sedflux_r;
	sedflux_r=(double*)calloc(number_of_gscl, sizeof(double));

	double *sedflux_l;
	sedflux_l=(double*)calloc(number_of_gscl, sizeof(double));

	flowpath::iterator i;
	int j=0;
	for (i = current_flowpath->begin(); i != coastline_cell_index; i++,j++)
	{
		curve_width = storage(i, discharge,j,sim_time);
	//	cout<<i->col<<" "<<"curve_width.right = "<<curve_width.right<<" "<<"curve_width.left = "<<curve_width.left<<endl;

		// define the erf of the actual flowpath cell DEBUGGEN RIEN, dit klopt niet 26 sept
		double cellpart_prev;
		double delta_cellpart = 2 * 1.0 / (double)(curve_width.right - i->col);
		double cellpart_current = delta_cellpart;

	//	double surf_facfluv_0 = erf (cellpart_current);
	//	cout<<j<<" "<<"deltacellpart "<< delta_cellpart<<" surf_facfluv_0= " <<surf_facfluv_0<<endl;
		
		// DIT IS DE AFZETTING IN HET FLOWPATH
	//	double sedflux_sum=0;
	//	for(int n=0;n<number_of_gscl;n++)
	//	{
	//		space[sim_time][i->row][i->col][n] +=sedflux_depo[j][n]*surf_facfluv_0*dt*dx;
	//		sedflux_sum +=sedflux_depo[j][n];
	//	}
	//		space[sim_time][i->row][i->col][0] += sedflux_sum*surf_facfluv_0;


		//righthand side		
				
		for(int l=i->col; l<curve_width.right;l++)
		{	

			//width is scaled to 0-2 domain appropriate for erf()
			cellpart_prev = cellpart_current;			
			cellpart_current += delta_cellpart;

			double surf_fac_fluv=0;
			surf_fac_fluv = (erf(cellpart_current) - erf(cellpart_prev));

			//cout<<"right "<<l<<" "<<"surf_fac_fluv = "<<surf_fac_fluv<<endl;
			
			double sedflux_sum_r=0;
			for(int n=0;n<number_of_gscl;n++)
			{
				sedflux_r[n] = 0.5 * sedflux_depo[j][n];
				
			//	sedflux_r[n] = ((double)(curve_width.right-i->col)/(curve_width.right-curve_width.left))* sedflux_depo[j][n];
			//	space[sim_time][i->row][l][n] +=sedflux_depo[j][n]*surf_fac_fluv*dt*dx;
			//	sedflux_sum_r +=sedflux_depo[j][n]*surf_fac_fluv*dt*dx;
				sedflux_sum_r +=sedflux_r[n]*surf_fac_fluv*dt*dx;
			}
				space[sim_time][i->row][l][0] += sedflux_sum_r;

			// to divide the sediment from coarse (graincell_cl =0) to fine (graincell_cl = num_of_gscl)
			// assumption made: within a single timestep the grainsizes-classes stack and don't mix

				int graincell_cl =0;
				double height = 0.0;
				while (height < sedflux_sum_r) {
					if (height + sedflux_r[graincell_cl] < sedflux_sum_r) {
						space[sim_time][i->row][l][graincell_cl+1]
							+= sedflux_r[graincell_cl];
						height += sedflux_r[graincell_cl++];
					} else {
						space[sim_time][i->row][l][graincell_cl+1]
							+= sedflux_sum_r - height;
						sedflux_r[graincell_cl] -= sedflux_sum_r - height;
						height = sedflux_sum_r;
					}

				}
		
		}//end rightside

		delta_cellpart = 2 * 1 / (double)(i->col - curve_width.left);
		cellpart_current = delta_cellpart;
		// lefthand side
		for(int m=(i->col)-1; m>curve_width.left; m--)
		{
			cellpart_prev = cellpart_current;			
			cellpart_current += delta_cellpart;
			double surf_fac_fluv = (erf(cellpart_current) - erf(cellpart_prev));
			
		//	cout<<"left "<<l<<" "<<"surf_fac_fluv = "<<surf_fac_fluv<<endl;
			double sedflux_sum_l=0;
			for(int n=0;n<number_of_gscl;n++)
			{				
				sedflux_l[n] = 0.5 *sedflux_depo[j][n];
			//	space[sim_time][i->row][m][n] +=sedflux_depo[j][n]*surf_fac_fluv*dt*dx;
			//	sedflux_sum_l +=sedflux_depo[j][n]*surf_fac_fluv*dx*dt;
				sedflux_sum_l +=sedflux_l[n]*surf_fac_fluv*dx*dt;
			}

			space[sim_time][i->row][m][0] += sedflux_sum_l*surf_fac_fluv;

			int graincell_cl=0;
			double height =0.0;
			// to divide the sediment from coarse (graincell_cl =0) to fine (graincell_cl = num_of_gscl)
			// assumption made: grainsizes-classes stack and don't mix
				
				while (height < sedflux_sum_l) {
					if (height + sedflux_l[graincell_cl] < sedflux_sum_l) 
					{
						space[sim_time][i->row][m][graincell_cl+1]
							+= sedflux_l[graincell_cl];
						height += sedflux_l[graincell_cl++];
					} 
					
					else 
					{
						space[sim_time][i->row][m][graincell_cl+1]
							+= sedflux_sum_l - height;
						sedflux_l[graincell_cl] -= sedflux_sum_l - height;
						height = sedflux_sum_l;
					}

				}
		}// end leftside

	}// end flowpath


}*/

void fluv_lat(flowpath* current_flowpath,flowpath::iterator coastline_cell_index, long sim_time, double *discharge, double **sedflux_depo)

{

	// to divide the depositional sediment flux over a left and righthand component
	double *sedflux_r;
	sedflux_r=(double*)calloc(number_of_gscl, sizeof(double));

	double *sedflux_l;
	sedflux_l=(double*)calloc(number_of_gscl, sizeof(double));

	double *l_scale;
	l_scale =(double *)calloc(number_of_colums, sizeof(double));
	
	flowpath::iterator i;
	int j=0;
	for (i = current_flowpath->begin(); i != coastline_cell_index; i++,j++)
	{
		int L_fluv = 8; //this determines the width of fluvial deposition zone
						// here it is done statically but should be dynamically in future

		int lat_width_r=0; // determines the actual sedimentation zone width (right)
		int lat_width_l=0; // determines the actual sedimentation zone width (left)

		
		//righthand side
		
		// boundary condition
		// latwidth should not exceed grid boundaries
			lat_width_r = i->col + L_fluv;
			if (lat_width_r > number_of_colums) {lat_width_r = number_of_colums;}

		double sum=0;

		int o =1;
		for(int l=i->col; l<lat_width_r;l++, o++)
		{
		
			l_scale[o]=(double)(l-i->col)/(L_fluv);
			double surf_fac_r = (erf(l_scale[o])-erf(l_scale[o-1])) ;
			sum=sum+ surf_fac_r;

			//cout<< "l = "<<l<<" surf_fac_r = "<<surf_fac_r<<" "<<sum<<endl;
			double sedflux_sum_r=0;
			for(int n=0;n<number_of_gscl;n++)
			{
				sedflux_r[n] = 0.5 * sedflux_depo[j][n];
			//	cout<<" n "<<n<<" sedflux_r_n "<<sedflux_r[n]<<endl;
				sedflux_sum_r = sedflux_sum_r+sedflux_r[n]*surf_fac_r*dt*dx;
			}
				space[sim_time][i->row][l][0] =space[sim_time][i->row][l][0]+ sedflux_sum_r;
			
			// to divide the sediment from coarse (graincell_cl =0) to fine (graincell_cl = num_of_gscl)
			// assumption made: within a single timestep the grainsizes-classes stack and don't mix

				int graincell_cl =0;
				double height = 0.0;
				while (height < sedflux_sum_r) 
				{
					if (height + sedflux_r[graincell_cl] < sedflux_sum_r) 
					{
						space[sim_time][i->row][l][graincell_cl+1]+=sedflux_r[graincell_cl];
							height+=sedflux_r[graincell_cl];
						graincell_cl++;
					} 
					else 
					{
						space[sim_time][i->row][l][graincell_cl+1]+=sedflux_sum_r - height;
						sedflux_r[graincell_cl] =sedflux_r[graincell_cl] -(sedflux_sum_r - height);
						height = sedflux_sum_r;
					}
				}//end while

			

		}// end right side


		// left hand side

		// boundary condition
		// latwidth should not exceed grid boundaries
		lat_width_r = i->col - L_fluv;
		if (lat_width_r < 1) {lat_width_r = 1;}

		int graincell_cl = 0;
		sum = 0;

		int p =1;
		for(int l=i->col+1; l>lat_width_r;l--, p++)
		{

			l_scale[p]=(double)(i->col-l)/(L_fluv);
			double surf_fac_l = (erf(l_scale[p])-erf(l_scale[p-1])) ;
			sum=sum+ surf_fac_l;

			double sedflux_sum_l=0;
			for(int n=0;n<number_of_gscl;n++)
			{
				sedflux_l[n] = 0.5 * sedflux_depo[j][n];
			//	cout<<" n "<<n<<" sedflux_r_n "<<sedflux_r[n]<<endl;
				sedflux_sum_l = sedflux_sum_l+sedflux_l[n]*surf_fac_l*dt*dx;
			}
				space[sim_time][i->row][l][0] =space[sim_time][i->row][l][0]+ sedflux_sum_l;
			
				// to divide the sediment from coarse (graincell_cl =0) to fine (graincell_cl = num_of_gscl)
			// assumption made: within a single timestep the grainsizes-classes stack and don't mix

				int graincell_cl =0;
				double height = 0.0;
				while (height < sedflux_sum_l) 
				{
					if (height + sedflux_l[graincell_cl] < sedflux_sum_l) 
					{
						space[sim_time][i->row][l][graincell_cl+1]+=sedflux_l[graincell_cl];
							height+=sedflux_l[graincell_cl];
						graincell_cl++;
					} 
					else 
					{
						space[sim_time][i->row][l][graincell_cl+1]+=sedflux_sum_l - height;
						sedflux_l[graincell_cl] =sedflux_l[graincell_cl] -(sedflux_sum_l - height);
						height = sedflux_sum_l;
					}
				}//end while







		}// end for left hand side
	}
}
