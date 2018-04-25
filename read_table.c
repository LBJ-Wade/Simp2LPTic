#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "allvars.h"
#include "proto.h"

/*! \file read_table.c
 *  \brief to get various kinds of tables calculated beforehand for different DEDMI models
 */
#ifdef HUBBLE_TABLE 
void me_init_hubble_table(void)
{
  FILE *fd;
  int loop;
  double a[HUBBLE_TABLE_LENGTH],h[HUBBLE_TABLE_LENGTH];
//  double drift[HUBBLE_TABLE_LENGTH],gravkick[HUBBLE_TABLE_LENGTH],hydrokick[HUBBLE_TABLE_LENGTH];
  char fname[500];
  strcpy(fname,FileWithInputHubble);
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
	printf("reading `%s' ...\n", fname);
    fflush(stdout);
  printf("loading Hubble Table ...\n");  
  for(loop=0;loop<HUBBLE_TABLE_LENGTH;loop++)
    {
		fscanf(fd,"%lf,%lf\n",&a[loop],&h[loop]);	
		//printf("input, a=%lf,h=%lf\n",a[loop],h[loop]); //debug use
//		drift[loop] = 1.0/(All.Hubble*h[loop]*a[loop]*a[loop]*a[loop]);
//		gravkick[loop] = 1.0/(All.Hubble*h[loop]*a[loop]*a[loop]);
//		hydrokick[loop] = 1.0/(All.Hubble*h[loop]*pow(a[loop], 3 * GAMMA_MINUS1) * a[loop]);
	}
  MeHubbleAcc = gsl_interp_accel_alloc();
  MeHubbleSpline = gsl_spline_alloc(gsl_interp_cspline, HUBBLE_TABLE_LENGTH);
//  MeHubbleSplineDRIFT = gsl_spline_alloc(gsl_interp_cspline, HUBBLE_TABLE_LENGTH);
//  MeHubbleSplineGRAVKICK = gsl_spline_alloc(gsl_interp_cspline, HUBBLE_TABLE_LENGTH);
//  MeHubbleSplineHYDROKICK = gsl_spline_alloc(gsl_interp_cspline, HUBBLE_TABLE_LENGTH);
  gsl_spline_init(MeHubbleSpline, a, h, HUBBLE_TABLE_LENGTH);
//  gsl_spline_init(MeHubbleSplineDRIFT, a, drift, HUBBLE_TABLE_LENGTH);
//  gsl_spline_init(MeHubbleSplineGRAVKICK, a, gravkick, HUBBLE_TABLE_LENGTH);
//  gsl_spline_init(MeHubbleSplineHYDROKICK, a, hydrokick, HUBBLE_TABLE_LENGTH);
  printf("Hubble Table loading done.\n");
  fflush(stdout);
}
#endif
#ifdef DMMASS_TABLE
void me_init_dmmass_table(void)
{
  FILE *fd;
  int loop;
  double a[DMMASS_TABLE_LENGTH],m[DMMASS_TABLE_LENGTH];
  char fname[100];
  strcpy(fname,"dmmass_table.txt");
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
	printf("reading `%s' ...\n", fname);
    fflush(stdout);
  printf("loading DMMass Table ...\n");  
  for(loop=0;loop<DMMASS_TABLE_LENGTH;loop++)
    {
		fscanf(fd,"%lf,%lf\n",&a[loop],&m[loop]);
		//printf("mass=%f\n",m[loop]);	
	}
  MeDMMassAcc = gsl_interp_accel_alloc();
  MeDMMassSpline = gsl_spline_alloc(gsl_interp_cspline, DMMASS_TABLE_LENGTH);
  gsl_spline_init(MeDMMassSpline, a, m, DMMASS_TABLE_LENGTH);
  printf("DMMass Table loading done.\n");
  fflush(stdout);
}
#endif
#ifdef DRAG_TABLE
void me_init_drag_table(void)
{
  FILE *fd;
  int loop;
  double a[DRAG_TABLE_LENGTH],d[DRAG_TABLE_LENGTH];
  char fname[100];
  strcpy(fname,"drag_table.txt");
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
	printf("reading `%s' ...\n", fname);
    fflush(stdout);
  printf("loading DRAG Table ...\n");  
  for(loop=0;loop<DRAG_TABLE_LENGTH;loop++)
    {
		fscanf(fd,"%lf,%lf\n",&a[loop],&d[loop]);
		//printf("drag=%f\n",d[loop]);	
	}
  All.MeDRAGAcc = gsl_interp_accel_alloc();
  All.MeDRAGSpline = gsl_spline_alloc(gsl_interp_cspline, DRAG_TABLE_LENGTH);
  gsl_spline_init(All.MeDRAGSpline, a, d, DRAG_TABLE_LENGTH);
  printf("DRAG Table loading done.\n");
  fflush(stdout);
}
#endif
#ifdef DE_TABLE
void me_init_de_table(void)
{
  FILE *fd;
  int i,j,dummy;
  double a[DE_TABLE_LENGTH_A],k[DE_TABLE_LENGTH_K];
  double *data = malloc(DE_TABLE_LENGTH_A * DE_TABLE_LENGTH_K * sizeof(double));
  double *r = malloc(DE_TABLE_LENGTH_A * DE_TABLE_LENGTH_K * sizeof(double));
  char fname[100]; 
  strcpy(fname,"de_table.txt");
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
	printf("reading `%s' ...\n", fname);
    fflush(stdout);
  printf("loading DE Table ...\n");  
  
  fscanf(fd,"%lf",&dummy);
  for(i=0;i<DE_TABLE_LENGTH_K;i++)
  {
	fscanf(fd,",%lf",&k[i]);
  }
  fscanf(fd,"\n");
  for(i=0;i<DE_TABLE_LENGTH_A;i++)
  {
	  fscanf(fd,"%lf",&a[i]);
	  for(j=0;j<DE_TABLE_LENGTH_K;j++)
	  {
		  fscanf(fd,",%lf",&data[i*DE_TABLE_LENGTH_K+j]);
	  }
	  fscanf(fd,"\n");
  }
  for(i=0;i<DE_TABLE_LENGTH_A;i++)
  {
	  for(j=0;j<DE_TABLE_LENGTH_K;j++)
	  {
		  gsl_spline2d_set(All.MeDESpline,r,a[i],k[j],data[i*DE_TABLE_LENGTH_K+j]);
	  }
  }
  gsl_spline2d_init(All.MeDESpline, a, k, r, DE_TABLE_LENGTH_A, DE_TABLE_LENGTH_K);
  All.MeDEAcca = gsl_interp_accel_alloc();
  All.MeDEAcck = gsl_interp_accel_alloc();
  All.MeDESplinea = gsl_spline_alloc(gsl_interp_cspline, DE_TABLE_LENGTH);
  All.MeDESplineb = gsl_spline_alloc(gsl_interp_cspline, DE_TABLE_LENGTH);
  printf("DE Table loading done.\n");
  fflush(stdout);
}
#endif
