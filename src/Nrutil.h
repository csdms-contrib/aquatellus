#ifndef _NRUTIL_H_
#define _NRUTIL_H_

#include <stdlib.h>
#include <stdio.h>

void
nrerror (const char error_text[])
//char error_text[];
{
  //void exit();

  fprintf (stderr, "Numerical Recipes run-time error...\n");
  fprintf (stderr, "%s\n", error_text);
  fprintf (stderr, "...now exiting to system...\n");
  //exit(1);
}

double *
vect (int nl, int nh)
//int nl,nh;
{
  double *v;

  v = (double *)malloc ((unsigned)(nh - nl + 1) * sizeof (double));
  if (!v)
    nrerror ("allocation failure in vect()");
  return v - nl;
}

int *
ivect (int nl, int nh)
//int nl,nh;
{
  int *v;

  v = (int *)malloc ((unsigned)(nh - nl + 1) * sizeof (int));
  if (!v)
    nrerror ("allocation failure in ivect()");
  return v - nl;
}

double *
dvect (int nl, int nh)
//int nl,nh;
{
  double *v;

  v = (double *)malloc ((unsigned)(nh - nl + 1) * sizeof (double));
  if (!v)
    nrerror ("allocation failure in dvect()");
  return v - nl;
}

double **
matrix (int nrl, int nrh, int ncl, int nch)
//int nrl,nrh,ncl,nch;
{
  int i;

  double **m;

  m = (double **)malloc ((unsigned)(nrh - nrl + 1) * sizeof (double *));
  if (!m)
    nrerror ("allocation failure 1 in matrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (double *)malloc ((unsigned)(nch - ncl + 1) * sizeof (double));
    if (!m[i])
      nrerror ("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  return m;
}

double **
dmatrix (int nrl, int nrh, int ncl, int nch)
//int nrl,nrh,ncl,nch;
{
  int i;

  double **m;

  m = (double **)malloc ((unsigned)(nrh - nrl + 1) * sizeof (double *));
  if (!m)
    nrerror ("allocation failure 1 in dmatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (double *)malloc ((unsigned)(nch - ncl + 1) * sizeof (double));
    if (!m[i])
      nrerror ("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}

int **
imatrix (int nrl, int nrh, int ncl, int nch)
//int nrl,nrh,ncl,nch;
{
  int i,
  **m;

  m = (int **)malloc ((unsigned)(nrh - nrl + 1) * sizeof (int *));
  if (!m)
    nrerror ("allocation failure 1 in imatrix()");
  m -= nrl;

  for (i = nrl; i <= nrh; i++)
  {
    m[i] = (int *)malloc ((unsigned)(nch - ncl + 1) * sizeof (int));
    if (!m[i])
      nrerror ("allocation failure 2 in imatrix()");
    m[i] -= ncl;
  }
  return m;
}

double **
submatrix (double **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl,
           int newcl)
//double **a;
//int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
  int i,
    j;

  double **m;

  m = (double **)malloc ((unsigned)(oldrh - oldrl + 1) * sizeof (double *));
  if (!m)
    nrerror ("allocation failure in submatrix()");
  m -= newrl;

  for (i = oldrl, j = newrl; i <= oldrh; i++, j++)
    m[j] = a[i] + oldcl - newcl;

  return m;
}

void
free_vect (double *v, int nl, int nh)
//double *v;
//int nl,nh;
{
  free ((char *)(v + nl));
}

void
free_ivect (int *v, int nl, int nh)
//int *v,nl,nh;
{
  free ((char *)(v + nl));
}

void
free_dvect (double *v, int nl, int nh)
//double *v;
//int nl,nh;
{
  free ((char *)(v + nl));
}

void
free_matrix (double **m, int nrl, int nrh, int ncl, int nch)
//double **m;
//int nrl,nrh,ncl,nch;
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((char *)(m[i] + ncl));
  free ((char *)(m + nrl));
}

void
free_dmatrix (double **m, int nrl, int nrh, int ncl, int nch)
//double **m;
//int nrl,nrh,ncl,nch;
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((char *)(m[i] + ncl));
  free ((char *)(m + nrl));
}

void
free_imatrix (int **m, int nrl, int nrh, int ncl, int nch)
//int **m;
//int nrl,nrh,ncl,nch;
{
  int i;

  for (i = nrh; i >= nrl; i--)
    free ((char *)(m[i] + ncl));
  free ((char *)(m + nrl));
}

void
free_submatrix (double **b, int nrl, int nrh, int ncl, int nch)
//double **b;
//int nrl,nrh,ncl,nch;
{
  free ((char *)(b + nrl));
}

double **
convert_matrix (double *a, int nrl, int nrh, int ncl, int nch)
//double *a;
//int nrl,nrh,ncl,nch;
{
  int i,
    j,
    nrow,
    ncol;

  double **m;

  nrow = nrh - nrl + 1;
  ncol = nch - ncl + 1;
  m = (double **)malloc ((unsigned)(nrow) * sizeof (double *));
  if (!m)
    nrerror ("allocation failure in convert_matrix()");
  m -= nrl;
  for (i = 0, j = nrl; i <= nrow - 1; i++, j++)
    m[j] = a + ncol * i - ncl;
  return m;
}

void
free_convert_matrix (double **b, int nrl, int nrh, int ncl, int nch)
//double **b;
//int nrl,nrh,ncl,nch;
{
  free ((char *)(b + nrl));
}

#endif /* _NRUTIL_H_ */
