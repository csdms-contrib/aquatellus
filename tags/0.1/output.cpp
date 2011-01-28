//      this function writes simple output to files
// 1.   a SURFER MAP of the topographical heights at a certain time step

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "main.h"
#include "Space.h"
#include "flowpath.h"
#include "inputcontrols.h"
#include "output.h"

using namespace std;


#define FILE_NAME1 "topo.grd"
#define FILE_NAME2 "flowpath.dat"
//int num_outp = 5; // this is the number of times output is wanted.



void
open_data_files (FILE ** ftopo, FILE ** fpath, int time)
{
  char filename[100];
  char *filePtr = &filename[0];
  ::sprintf (filePtr, "topo%d.grd", time);
  if ((*ftopo = fopen (filename, "w")) == NULL)
  {
    cerr << "error while opening file" << endl;
    exit (1);
  }

  ::sprintf (filePtr, "flowpath%d.dat", time);
  if ((*fpath = fopen (filename, "w")) == NULL)
  {
    cerr << "error while opening file" << endl;
    exit (1);
  }
}

void
close_data_files (FILE ** ftopo, FILE ** fpath)
{
  fclose (*ftopo);
  fclose (*fpath);
}

// 1. here the header for the SURFER GRID FILE is written
void
write_topographical_map_header (FILE ** ftopo)
{
  // to indicate the data format as ASCII
  fprintf (*ftopo, "DSAA \n");
  // the total number of colums (x-axis) and rows(y-axis) are defined for drawing
  fprintf (*ftopo, "%3d %3d\n", number_of_colums, number_of_rows);
  // here in the header the minimum and maximum x-values are defined for drawing
  fprintf (*ftopo, "  0  %5d \n", number_of_colums);
  // here in the header the minimum and maximum y-values are defined for drawing
  fprintf (*ftopo, "  0  %5d \n", number_of_rows);
  // here in the header the minimum and maximum z-values are defined for drawing
  fprintf (*ftopo, " 40  110 \n");
}


void
write_topographical_map (FILE ** ftopo, long sim_time)
{
  int i, j;
  for (i = 0; i < number_of_rows; i++)
  {
    for (j = 0; j < number_of_colums; j++)
    {
      double height = 0;
      for (int t = 0; t <= sim_time; t++)
      {
        height += space[t][i][j][0];

      }
      fprintf (*ftopo, "%5.5f ", height);
    }

    fprintf (*ftopo, "\n");
  }
}


// temporary calculation of weighted average of grainsize distribution
double **
calculate_grain_distr (double *grain_size, int Xpos)
{
  int i, j, k;
  double **median;
  int nx = number_of_colums;
  // code has to be totally adapted for longitudinal Xsection
  // int nx = number_of_rows;// in case of longitudinal crosssection

  // allocation of median
  if ((median =
       (double **) calloc ((end_of_times / dt), sizeof (double *))) == NULL)
    printf ("median array could not be allocated");
  for (i = 0; i < (end_of_times / dt); i++)
  {
    if ((median[i] = (double *) calloc (nx, sizeof (double))) == NULL)
      printf ("median array could not be allocated");
  }


  double tmpsum;
  for (j = 0; j < end_of_times / dt; j++)       //Calculate weighted average grainsize using the farmer-woodenshoes method
  {
    for (i = 0; i < nx; i++)
    {
      if (space[j][Xpos][i][0] != 0)
      {
        tmpsum = 0;
        for (k = 1; k <= number_of_gscl; k++)
        {
          tmpsum += (space[j][Xpos][i][k] * grain_size[k - 1]);
          median[j][i] = 10000 * (tmpsum);      // in micron


        }

      }
      else
        median[j][i] = 4;
      //cout<<" time " <<j<<" position "<<i<<" "<<median[j][i]<<endl; 
    }
  }
  return median;
}


// this function is necessary for the well routine
// temporary calculation of weighted average of grainsize distribution
double *
calculate_grain_distr_well (double *grain_size, int wellpos_row,
                            int wellpos_col)
{
  int j, k;
  double *median_point;

  // allocation of median
  if ((median_point =
       (double *) calloc ((end_of_times / dt), sizeof (double))) == NULL)
    printf ("median array could not be allocated");

  double tmpsum;
  for (j = 0; j < end_of_times / dt; j++)       //Calculate weighted average grainsize using the farmer-woodenshoes method
  {
    if (space[j][wellpos_row][wellpos_col][0] != 0)
    {
      tmpsum = 0;
      for (k = 1; k <= number_of_gscl; k++)
      {
        tmpsum +=
          (space[j][wellpos_row][wellpos_col][k - 1] * grain_size[k - 1]);
        median_point[j] = 10000 * (tmpsum);     ///space[j][wellpos_row][wellpos_col][0]);  // in micron
        cout << "time " << j << " " << tmpsum << " depth " <<
          space[j][wellpos_row][wellpos_col][0] << " median " <<
          median_point[j] << endl;
      }

    }
    else
      median_point[j] = 4;
    //cout<<" time " <<j<<" "<<median_point[j]<<endl;       
  }

  return median_point;
}

void
write_median_psfile (double **median, int Xpos)
{

  long i, j;
  short Lmarg, Rmarg, Tmarg, Bmarg;     // left, right, top and bottem margins
  long hhmin, hhmax;
  double hmin, hmax, scx, scy, huemedian, medianmin, medianmax, medianminn,
    medianmaxx, dmedian;
  FILE *psfp;
  char f_name[255];

  // these determine the longitudinal scale for output
  //x-pos start &stop ps-output
  short xstrtps = 10;
  short xstopps = number_of_colums - 1;

  Rmarg = 10;
  Tmarg = 120;
  Bmarg = 120;                  /* 0.5 inch margins */
  Lmarg = 30;                   /* 30 pts for vertical scale bar        */

  strcpy (f_name, "xsection.ps");

  psfp = fopen (f_name, "w");

  fprintf (psfp, "%%!\n");
  fprintf (psfp, "gsave\n");
  fprintf (psfp, "/Times-Bold findfont 10 scalefont setfont\n");
  fprintf (psfp, "grestore\n");

  /* draw grainsize color coding bar */
  /* first determine grainsize range */
  medianmax = -1000000;
  medianmin = 1000000;

  for (i = 0; i < (end_of_times / dt); i++)
  {
    for (j = xstrtps; j < xstopps; j++)
    {
      medianmin = 0;            //medianminn = MIN(median[i][j],medianmin);
      medianmax = 400;          //medianmaxx = MAX(median[i][j],medianmax); 
    }
  }


  dmedian = (medianmax - medianmin) / 10;

  fprintf (psfp, "gsave\n");
  fprintf (psfp, "%d %d translate\n", Lmarg, Bmarg);
  fprintf (psfp, "/Helvetica findfont 10 scalefont setfont\n");
  fprintf (psfp,
           "100 -8 moveto 0 0 0 setrgbcolor (Median Grainsize [mu]) show\n");
  fprintf (psfp, "grestore\n");

  fprintf (psfp, "gsave\n");
  fprintf (psfp, "%d %d translate\n", Lmarg + 25, Bmarg + 25);
  fprintf (psfp, "/Helvetica findfont 10 scalefont setfont\n");

  // Now draw 10 median classes - not log median  
  fprintf (psfp, "0 0 0 setrgbcolor\n");
  fprintf (psfp, "10 15 moveto\n0 -20 rlineto\nstroke\n");      //the first 

  for (i = 0; i <= 10; i++)
  {
    huemedian = (0.70 - 0.70 * i / 9);  // hue runs from 0(red) to 0.75(blue)
    fprintf (psfp, "%d 0 moveto\n", (i * 25) + 10);
    fprintf (psfp, "0 15 rlineto\n25 0 rlineto\n0 -15 rlineto\n-25 0 rlineto\n");       //boxfilldivider
    fprintf (psfp, "%5.3f 1.0 1.0 sethsbcolor\n fill\n", huemedian);    //COLOUR
    fprintf (psfp, "0 0 0 setrgbcolor\n");
    fprintf (psfp, "%d 0 moveto\n", (i * 25) + 10);
    fprintf (psfp, "0 15 rlineto\n25 0 rlineto\n0 -15 rlineto\n-25 0 rlineto\nstroke\n");       //box
    fprintf (psfp, "%d 0 moveto\n0 -5 rlineto\nstroke\n", ((i + 1) * 25) + 10); //divider line
    fprintf (psfp, "%d 0 moveto\n-10 -15 rmoveto\n (%8.0f) show\n", (i * 25) + 10, medianmin + i * dmedian);    //value

  }

  fprintf (psfp, "grestore\n");
  fprintf (psfp, "0.15 setlinewidth\n");
  fprintf (psfp, "0.0 setgray\n");

  /* Determine x & y scale */
  hmax = hhmax = -1000000;
  hmin = hhmin = 1000000;
  for (i = 0; i < (end_of_times / dt); i++)
  {
    for (j = xstrtps; j < xstopps; j++)
    {
      hmin = 82;                //MIN(get_topo_height(i, Xpos, j), hmin);
      hmax = 92;                //MAX(get_topo_height(i, Xpos, j),hmax);
    }
  }

  //cout<<"hmin and hmax "<<hmin<< " "<<hmax<<endl;

  scy = (612.0 - 60 - (Tmarg + Bmarg)) / (hmax - hmin);
  scx = (792.0 - (Lmarg + Rmarg)) / (xstopps - xstrtps);

  //this is for a section perpendicular to flowpath
  // draw vertical scale bar
  hhmin = hmin;
  hhmax = hmax;
  fprintf (psfp, "gsave\n");
  fprintf (psfp, "%d %d translate\n", -20, 0);  /* draw scale bar left of Lmarg */

  fprintf (psfp, "%d %d moveto\n", Lmarg, Bmarg);
  fprintf (psfp, "%d %3.0f lineto\n", Lmarg, (Bmarg + scy * (hhmax - hhmin)));
  fprintf (psfp, "stroke\n");
  fprintf (psfp, "/Helvetica findfont 6 scalefont setfont\n");

  for (i = hhmin; i <= hhmax; i++)
  {
    fprintf (psfp, "%d %3.0f moveto\n", Lmarg, Bmarg + scy * (i - hhmin));
    if (i % 1 == 0)
    {
      fprintf (psfp, "10 0 rlineto\n");
      fprintf (psfp, "stroke\n");
      fprintf (psfp, "%d %10.0f moveto\n", Lmarg, Bmarg + scy * (i - hhmin));
      fprintf (psfp, "13 -3 rmoveto\n");
      fprintf (psfp, "(%d) show\n", i);
    }
    else
    {
      fprintf (psfp, "0.15 setlinewidth\n");
      fprintf (psfp, "0.0 setgray\n");
      fprintf (psfp, "5 0 rlineto\n");
      fprintf (psfp, "stroke\n");
    }
  }
  fprintf (psfp, "grestore\n");

  for (i = 0; i < (end_of_times / dt); i++)
  {
    for (j = xstrtps; j < xstopps; j++)
    {

      huemedian = 0.70 - 0.70 * (median[i][j] - medianmin) / (medianmax - medianmin);   /* COLOUR 0.7=rood ..0.35=groen .. 0.0=blauw */

      if (i > 0)
      {
        /*      Draw a cell. nodestrat[i][j] is upper left corner,nodestrat[i+1][j+1] the lower right. */
        if (space[i][Xpos][j][0] < EPSloc)      /* zero thickness -> borrow median from neighbour */
          huemedian = 0.70 - 0.70 * (median[i][j + 1] - medianmin) / (medianmax - medianmin);   //COLOUR

        /*hier aangepast omdat topografische hoogte gebruikt moet worden?
           fprintf(psfp,"%5.3f 1.0 1.0 sethsbcolor\n",huemedian); //COLOUR
           fprintf(psfp,"%5.1f %5.1f moveto\n",Lmarg+scx*(j-xstrtps),  Bmarg+scy*(space[i][Xpos][j][0]-hmin));
           fprintf(psfp,"%5.1f %5.1f lineto\n",Lmarg+scx*(j+1-xstrtps),Bmarg+scy*(space[i][Xpos][j+1][0]-hmin));
           fprintf(psfp,"%5.1f %5.1f lineto\n",Lmarg+scx*(j+1-xstrtps),Bmarg+scy*(space[i-1][Xpos][j+1][0]-hmin));
           fprintf(psfp,"%5.1f %5.1f lineto\n",Lmarg+scx*(j-xstrtps),  Bmarg+scy*(space[i-1][Xpos][j][0]-hmin));
           fprintf(psfp,"%5.1f %5.1f lineto\n",Lmarg+scx*(j-xstrtps),  Bmarg+scy*(space[i][Xpos][j][0]-hmin));  
           fprintf(psfp,"fill\n"); */

        //hier aangepast omdat topografische hoogte gebruikt moet worden?
        fprintf (psfp, "%5.3f 1.0 1.0 sethsbcolor\n", huemedian);       //COLOUR
        fprintf (psfp, "%5.1f %5.1f moveto\n", Lmarg + scx * (j - xstrtps),
                 Bmarg + scy * (get_topo_height (i, Xpos, j) - hmin));
        fprintf (psfp, "%5.1f %5.1f lineto\n",
                 Lmarg + scx * (j + 1 - xstrtps),
                 Bmarg + scy * (get_topo_height (i, Xpos, j + 1) - hmin));
        fprintf (psfp, "%5.1f %5.1f lineto\n",
                 Lmarg + scx * (j + 1 - xstrtps),
                 Bmarg + scy * (get_topo_height (i - 1, Xpos, j + 1) - hmin));
        fprintf (psfp, "%5.1f %5.1f lineto\n", Lmarg + scx * (j - xstrtps),
                 Bmarg + scy * (get_topo_height (i - 1, Xpos, j) - hmin));
        fprintf (psfp, "%5.1f %5.1f lineto\n", Lmarg + scx * (j - xstrtps),
                 Bmarg + scy * (get_topo_height (i, Xpos, j) - hmin));
        fprintf (psfp, "fill\n");

      }


      if (i == 0)               /* put tickmarks every j and labels every 5 j   */
      {
        fprintf (psfp, "gsave\n");
        fprintf (psfp, "0 0 0 setrgbcolor\n");
        fprintf (psfp, "0.15 setlinewidth\n");
        // again some expiri's with get topoheight
        fprintf (psfp, "%5.1f %5.1f moveto\n", Lmarg + scx * (j - xstrtps),
                 Bmarg + scy * (space[i][Xpos][j][0] - hmin));
        //fprintf(psfp,"%5.1f %5.1f moveto\n",Lmarg+scx*(j-xstrtps),get_topo_height(0,Xpos,j)-hmin);

        if (j % 1 == 0)
        {                       /*outline tickmarks according to number */
          fprintf (psfp, "0 -3 rlineto\n");
          fprintf (psfp, "-3 -10 rmoveto\n");

          if (j % 5 == 0)
          {
            fprintf (psfp, "/Helvetica findfont 6 scalefont setfont\n");
            fprintf (psfp, "(%d) show\n", j);
          }
        }

        fprintf (psfp, "stroke\n");
        fprintf (psfp, "grestore\n");
      }

    }
  }


  fprintf (psfp, "showpage\n");


  fclose (psfp);
}


void
write_simulatedwell_psfile (double *median_point, int wellpos_row,
                            int wellpos_col)
{
  short i;
  short Lmarg, Rmarg, Tmarg, Bmarg;
  double hhmin, hhmax;
  double hmin, hmax, scx, scy, Mmin, Mmax, Mminn, Mmaxx;
  double red, green, blue;
  FILE *psfp;
  char f_name[255];

  Rmarg = 0;
  Tmarg = 0;
  Bmarg = 0;
  Lmarg = 10;

  // 30 pts for vertical scale bar        
  strcpy (f_name, "well.ps");
  psfp = fopen (f_name, "w");

  fprintf (psfp, "%%!\n");
  fprintf (psfp, "gsave\n");

  // draw permeability color coding bar 
  // first determine perm range 
  Mmax = -1000000;
  Mmin = 1000000;

  for (i = 1; i < end_of_times / dt; i++)
  {
    Mmin = Mminn = MIN (median_point[i], Mmin);
    Mmax = Mmaxx = MAX (median_point[i], Mmax);
  }
  cout << " Mmax " << Mmax << endl;
  // Determine x & y scale 

  hmax = hhmax = -1000000;
  hmin = hhmin = 1000000;
  for (i = 1; i < end_of_times / dt; i++)
  {
    hmin = MIN (get_topo_height (i, wellpos_row, wellpos_col), hmin);
    hmax = MAX (get_topo_height (i, wellpos_row, wellpos_col), hmax);
  }
  cout << "hmin and hmax " << hmax - hmin << endl;

  scy = (480 - (Tmarg + Bmarg)) / (hmax - hmin);
  scx = (Lmarg + Rmarg);


  // Draw well FILL 
  fprintf (psfp, "0 0 0 setrgbcolor\n");
  fprintf (psfp, "0.02 setlinewidth\n");

  fprintf (psfp, "%5.1f %5.1f moveto\n", Lmarg + scx,
           Bmarg +
           scy *
           (get_topo_height (end_of_times / dt - 1, wellpos_row, wellpos_col)
            - hmin));
  fprintf (psfp, "%5.1f %5.1f lineto\n", Lmarg + scx,
           Bmarg + scy * (get_topo_height (0, wellpos_row, wellpos_col) -
                          hmin));
  fprintf (psfp, "stroke\n");   // the above 3 line plot the well position and the lowerboundary -'vertical'line
  fprintf (psfp, "%5.1f %5.1f moveto\n", Lmarg + scx, Bmarg + scy * (get_topo_height (0, wellpos_row, wellpos_col) - hmin));    //goto lower left position were the 'L'-line that has just been drawn

  for (i = 1; i < end_of_times / dt; i++)
  {
    if (space[i][wellpos_row][wellpos_col][0] != 0)     //to prevent it drawing zerothickness beds
    {
      red = 1;
      green = 1;
      blue = 0;

      fprintf (psfp, "%5.1f %5.1f lineto\n",
               (0.30 * (median_point[i])) + Lmarg + scx,
               Bmarg +
               scy * (get_topo_height (i - 1, wellpos_row, wellpos_col) -
                      hmin));
      fprintf (psfp, "%5.1f %5.1f lineto\n",
               (0.30 * (median_point[i])) + Lmarg + scx,
               Bmarg + scy * (get_topo_height (i, wellpos_row, wellpos_col) -
                              hmin));
      fprintf (psfp, "%5.1f %5.1f lineto\n", Lmarg + scx,
               Bmarg + scy * (get_topo_height (i, wellpos_row, wellpos_col) -
                              hmin));
      fprintf (psfp, "%5.1f %5.1f lineto\n", Lmarg + scx,
               Bmarg +
               scy * (get_topo_height (i - 1, wellpos_row, wellpos_col) -
                      hmin));
      fprintf (psfp, "%5.3f %5.3f %5.3f setrgbcolor\n fill\n", red, green,
               blue);
      fprintf (psfp, "stroke\n");

      fprintf (psfp, "0 0 0 setrgbcolor\n");
      fprintf (psfp, "%5.1f %5.1f moveto\n", Lmarg + scx,
               Bmarg +
               scy * (get_topo_height (i - 1, wellpos_row, wellpos_col) -
                      hmin));
      fprintf (psfp, "%5.1f %5.1f lineto\n",
               (0.30 * (median_point[i])) + Lmarg + scx,
               Bmarg +
               scy * (get_topo_height (i - 1, wellpos_row, wellpos_col) -
                      hmin));
      fprintf (psfp, "%5.1f %5.1f lineto\n",
               (0.30 * (median_point[i])) + Lmarg + scx,
               Bmarg + scy * (get_topo_height (i, wellpos_row, wellpos_col) -
                              hmin));
      fprintf (psfp, "%5.1f %5.1f lineto\n", Lmarg + scx,
               Bmarg + scy * (get_topo_height (i, wellpos_row, wellpos_col) -
                              hmin));
      fprintf (psfp, "%5.1f %5.1f lineto\n", Lmarg + scx,
               Bmarg +
               scy * (get_topo_height (i - 1, wellpos_row, wellpos_col) -
                      hmin));
      fprintf (psfp, "stroke\n");

      fprintf (psfp, "%5.1f %5.1f moveto\n", Lmarg + scx,
               Bmarg + scy * (get_topo_height (i, wellpos_row, wellpos_col) -
                              hmin));

    }
//              }// endif
//              fprintf(psfp,"stroke\n");
  }
  fprintf (psfp, "showpage\n");

  fclose (psfp);
}
