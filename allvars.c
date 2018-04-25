#include "allvars.h"
#include <gsl/gsl_spline.h>

#ifdef HUBBLE_TABLE
gsl_interp_accel *MeHubbleAcc;
gsl_spline *MeHubbleSpline;
#endif
#ifdef DMMASS_TABLE
gsl_interp_accel *MeDMMassAcc;
gsl_spline *MeDMMassSpline;
#endif
double HubbleUserA;
double OmegaUserA;

struct io_header_1 header1, header;

int WhichSpectrum;


int SphereMode;
int *Local_nx_table;

FILE *FdTmp, *FdTmpInput;

int Nmesh, Nsample;

long long IDStart;



char GlassFile[500];
char FileWithInputSpectrum[500];
#ifdef HUBBLE_TABLE
char FileWithInputHubble[500];
#endif
#ifdef DMMASS_TABLE
char FileWithInputDMmass[500];
#endif
int GlassTileFac;

double Box;
int Seed;

long long TotNumPart;

int NumPart;

int *Slab_to_task;

int NTaskWithN;

struct part_data *P;

int Nglass;

double InitTime;
double Redshift;
double MassTable[6];


char OutputDir[100], FileBase[100];
int NumFilesWrittenInParallel;


int ThisTask, NTask;

int Local_nx, Local_x_start;

int IdStart;

rfftwnd_mpi_plan Inverse_plan;
rfftwnd_mpi_plan Forward_plan;
unsigned int TotalSizePlusAdditional;
//fftw_real *Disp;
fftw_real *Workspace;
//fftw_complex *Cdata;


double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
double InputSpectrum_UnitLength_in_cm;
double G, Hubble;
double RhoCrit;

double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
double OmegaBaryon, HubbleParam;
double ShapeGamma;
double PrimordialIndex;
double Dplus;			/* growth factor */

#ifdef DIFFERENT_TRANSFER_FUNC
int Type, MinType, MaxType;
#endif

