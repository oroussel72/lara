#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
using namespace std;

#include "Global.h"
#include "Vector.h"
#include "Timer.h"
#include "Parameters.h"

/*
_______________________________________________________________________________

Default values
_______________________________________________________________________________

Geometry
_______________________________________________________________________________

*/


int Dimension = 1;
int ScaleNb = 0;
real XMin = -1.;
real XMax = 1.;
real YMin = -1.;
real YMax = 1.;
real ZMin = 0.;
real ZMax = 1.;

int NX = 0;
int NY = 1;
int NZ = 1;


int BCXMin = Periodic;
int BCXMax = Periodic; 
int BCYMin = Periodic;
int BCYMax = Periodic; 
int BCZMin = FreeSlip;
int BCZMax = FreeSlip;

/*
_______________________________________________________________________________

Multiresolution
_______________________________________________________________________________

*/

bool Multiresolution = false;
real Tolerance = 1.E-02;
int ScaleMin = 4;
bool UseVariableTolerance=false;
real MultiresolutionError = 5.;
bool RefineInCenter = false;
int RefineWidth = 20;

/*
_______________________________________________________________________________

Numerics
_______________________________________________________________________________

*/


int Equation = ConvectionDiffusion;
int Regime = Supersonic;
int Scheme = Roe;
int Reconstruction = MUSCL;
int Limiter = MinMod;
int Centering = ThirdOrder;
int WENOVariables = Characteristic;
int SpeciesNb = 1;
int StepNb = 2;
real TimeStep = 0.;
real CFL = 0.5;
real TotalTime = 2.;
int WriteNb = 20;
int WriteFormat = VTK;
int WriteType = ASCII;
bool WriteStrip = false;
bool WriteStripCoarse = false;
bool Restart = false;

/*
_______________________________________________________________________________

Physics
_______________________________________________________________________________

*/

bool IsInviscid = true;
real Celerity = 1.;
real Diffusivity = 0.;
real Viscosity = 1.716E-05;
bool ConstantViscosity = false;
Vector MassDiffusivity(1);
real Gamma = 1.4;
real R = 8.314;
real MolarMass = 2.9E-02;
real Pr = 0.71;
real Tr = 273.15;
real Ts = 110.4;

/*
_______________________________________________________________________________

Detonation
_______________________________________________________________________________

*/

bool Detonation = false;
real Ta = 25.;
real Ti = 0.3;
real Q0 = 25.;
real K0 = 164180.;
real RKFTolerance = 1.E-6;
real RKFTimeStepFactor = 1.E-3;

/*
_______________________________________________________________________________

Penalization
_______________________________________________________________________________

*/

bool UsePenalization = false;
real Permeability = 1.E-12;
real PenalTemperature = 273.15;
real PenalVelocity = 0.;
real Porosity = 1.;

/*
_______________________________________________________________________________

Internal parameters
_______________________________________________________________________________

*/

int RefreshNb = 0;
int IterationNo = 0;
int StepNo = 0;
real RemainingTime = 0.;
real AverageTimeStep = 0.;
int WriteNo = 0;
int QuantityNb = 1;
int StorageNb = 2;
bool Verbose = false;
real SpaceStep = 0.;
real Conductivity = 0.;
Vector MaxAverage;
int NodeNb = 0;
int LeafNb = 0;
real AverageNodeNb = 0.;
real AverageLeafNb = 0.;
real StartTime = 0.;

Vector PredictionError;
Vector PredictionErrorOld;
real PredictionFactor = 1.;
real InitialTolerance = 1.E-2;

Timer LaraTimer;

/*
_______________________________________________________________________________

Integral values
_______________________________________________________________________________

*/

real MaxEigenvalue=0.;
real KineticEnergy=0.;
real MaxDiffusivity=0.;
real TimeAverageMaxDensity=0.;


/*
_______________________________________________________________________________

Init parameters
_______________________________________________________________________________

*/

void initParameters()
{
  
#include "parameters.dat"  
  
  // Init RemainingTime
  
  RemainingTime = TotalTime;
  
  // If use penalization, use last species as mask
  
  if (UsePenalization)
    SpeciesNb++;
  
  // Compute NX, NY, NZ
  
  if (ScaleNb !=0)
  {
    NX = 1<<ScaleNb;
    NY = 1<<ScaleNb;
    NZ = 1<<ScaleNb;   
  }
   
  // Compute number of conservative quantities
  
  switch(Equation)
  {
    case ConvectionDiffusion:
      QuantityNb = 1;
      Celerity /= sqrt(Dimension);
      break;
      
    case Burgers:
      QuantityNb = Dimension;
      break;
      
    case NavierStokes:
      QuantityNb = Dimension + SpeciesNb + 1;
      break;
  };  

  // When an outlet boundary condition exists, impose StorageNb = 4
  
  if (BCXMin == Outlet || BCXMax == Outlet) StorageNb = 4;
  if (Dimension  > 1)
    if (BCYMin == Outlet || BCYMax == Outlet) StorageNb = 4;
  if (Dimension  > 2)
    if (BCZMin == Outlet || BCZMax == Outlet) StorageNb = 4;

  // For McCormack schemes, impose two steps
  
  if (Scheme == McCormack)
    StepNb = 2;

/*  
  // Compute StorageNb in function of StepNb
  
  if (StepNb < 4)
    StorageNb = 2;
  else
    StorageNb = 4;
*/

  // Compute conductivity parameters
  
  Conductivity = Gamma*R/(Pr*MolarMass*(Gamma-1.));
  
  // Adapt different parameters in function of QuantityNb
  
  MaxAverage.raise(QuantityNb);
  PredictionError.raise(QuantityNb);
  PredictionErrorOld.raise(QuantityNb);
  
  if (UsePenalization) 
  {
    MassDiffusivity.resize(SpeciesNb-1);
    MassDiffusivity.set(SpeciesNb-1,0.);
  }

  // Compute RefreshNb (in case it is 0)
  
  if (RefreshNb == 0)
  {
    if (Dimension == 1)
      RefreshNb = 10;
    else if (Dimension == 2)
      RefreshNb = 2;
    else
      RefreshNb = 1;
  }

  // Compute initial tolerance
  
  InitialTolerance = Tolerance;
}

/*
_______________________________________________________________________________

Initial condition
_______________________________________________________________________________

*/

Vector initCondition(const real x, const real y, const real z)
{
  Vector result(QuantityNb), pen(QuantityNb);
  real u, v, w;
  real rho, p, rhoE, T;
  Vector Y(SpeciesNb);
  int i;
  
  u=v=w=rho=p=rhoE=T=0.;
  
  // Include file containing the definition of the initial condition
  
#include "initial.dat"
 
  // For Navier-Stokes, compute rhoE and set the vector of conservative quantities
  
  switch(Equation)
  {
    case ConvectionDiffusion:
    default:
      result.set(1,u);
      break;
    
    case Burgers:
      result.set(1,u);
      if (Dimension > 1) result.set(2,v);
      if (Dimension > 2) result.set(3,w);
      break;
    
    case NavierStokes:
  
      if (rhoE == 0. && p != 0.)  // Initial condition with pressure
	rhoE = p/(Gamma-1.) + 0.5*rho*(u*u+v*v+w*w) + ((Detonation)? Q0*rho*Y.at(1):0.);
    
      if (rhoE == 0. && p == 0.) // Initial condition with temperature
      {
	p = rho*R*T/MolarMass;
	rhoE = p/(Gamma-1.) + 0.5*rho*(u*u+v*v+w*w) + ((Detonation)? Q0*rho*Y.at(1):0.);
      }
    
      result.set(1,rho);
      result.set(2,rho*u);
        
      if (Dimension > 1) result.set(3,rho*v);
      if (Dimension > 2) result.set(4,rho*w);
    
      result.set(Dimension+2,rhoE);
      
      if (SpeciesNb > 1)
      {
	for (i=1; i <= SpeciesNb-1; i++)
	  result.set(Dimension+2+i,rho*Y.at(i));
      }
      
      if (UsePenalization)
	result.set(QuantityNb, (isImmersed(x,y,z,0.)? 0.:1.));	
      
      break;      
  };
  
  return result;
}
  
/*
_______________________________________________________________________________

Zones to immerse
_______________________________________________________________________________

*/

bool isImmersed(const real x, const real y, const real z, const real t)
{
  bool immersed=false;
  
#include "immersed.dat"
  
  return immersed;
}


/*
_______________________________________________________________________________

Indicate time to show
_______________________________________________________________________________

*/

bool timeToShow()
{
  return ((IterationNo-1)%RefreshNb == 0);
}

/*
_______________________________________________________________________________

Show performance
_______________________________________________________________________________

*/

void showPerformance()
{
  int day, hour, min, sec;
  unsigned int rest; 
  double CPUtime, estimatedCPUtime, remainingCPUtime;
  FILE* f;
  int N = 1<<(ScaleNb*Dimension);
  
  if (ScaleNb == 0)
    N = NX*NY*NZ;
  
  // --- Compute CPU time
  
  CPUtime = LaraTimer.CPUTime();
  estimatedCPUtime = CPUtime/(TotalTime-StartTime-RemainingTime)*TotalTime;
  remainingCPUtime = CPUtime/(TotalTime-StartTime-RemainingTime)*RemainingTime;
  
  // Write on file
  
  f = fopen("performance.out","w");
  
  fprintf(f,"%23s : %15s\n", "Solver", (Multiresolution)? "MR":"FV");
  fprintf(f,"%23s : %15i\n", "Dimension", Dimension);

  if (ScaleNb != 0)
    fprintf(f,"%23s : %15i\n", "Scales", ScaleNb);
  
  fprintf(f,"%23s : %15i\n", "Cells", N);
  fprintf(f, "%23s : %15i\n", "Time accuracy order", StepNb);

  if (Multiresolution)
  {    
    fprintf(f, "%23s : %13.1f %1s\n", "Cell compression", (100.*AverageNodeNb)/N, "%");
    fprintf(f, "%23s : %13.1f %1s\n", "Leaf compression", (100.*AverageLeafNb)/N, "%");
  }  
  if (RemainingTime > 0.)
  {
    fprintf(f,"%23s : %13.1f %1s\n", "Percentage done", (TotalTime-RemainingTime)*100./TotalTime, "%");
    fprintf(f,"%23s : %15i\n", "Iterations (elapsed)", IterationNo);
    fprintf(f,"%23s : %15i\n", "Iterations (estimated)", (int)(IterationNo*TotalTime/(TotalTime-RemainingTime)));
    fprintf(f,"%23s : %15.8e\n", "Time step", TimeStep);
    fprintf(f,"%23s : %15.8e\n", "Physical time (elapsed)", TotalTime-RemainingTime);
  }
  else  
  {
    fprintf(f,"%23s : %15i\n", "Iterations", IterationNo);
    fprintf(f,"%23s : %15.8e\n", "Average time step", AverageTimeStep);
  }

  fprintf(f, "%23s : %15.8e\n", "Physical time (total)", TotalTime);
  fprintf(f, "%23s : %15.8e\n", "CPU time (s)", CPUtime);
  fprintf(f, "%23s : %15.8e\n\n", "CPU time / it. x pt", CPUtime/((1.*N)*IterationNo));
  
  // Write estimated CPU time 
    
  rest =  (unsigned int)(estimatedCPUtime);
  day  =  rest/86400;
  rest %= 86400;
  hour =  rest/3600;
  rest %= 3600;
  min	= rest/60;
  rest %= 60;
  sec  =  rest;
  rest =  (unsigned int)(estimatedCPUtime);
  
  if (rest >= 86400)
    fprintf(f, "%23s : %5d day %2d h %2d min %2d s\n", "Total CPU time", day, hour, min, sec);

  if ((rest < 86400)&&(rest >= 3600))
    fprintf(f, "%23s : %2d h %2d min %2d s\n", "Total CPU time", hour, min, sec);

  if ((rest < 3600)&&(rest >= 60))
    fprintf(f, "%23s : %2d min %2d s\n", "Total CPU time", min, sec);

  if (rest < 60)
    fprintf(f, "%23s : %2d s\n", "Total CPU time", sec);

  if (RemainingTime > 0.)
  {
    // Computing remaining time 
      
    rest =  (unsigned int)(remainingCPUtime);
    day  =  rest/86400;
    rest %= 86400;
    hour =  rest/3600;
    rest %= 3600;
    min	= rest/60;
    rest %= 60;
    sec  =  rest;
    rest =  (unsigned int)(remainingCPUtime);
   
    // Write remaining CPU time 
    if (rest >= 86400)
      fprintf(f, "%23s : %5d day %2d h %2d min %2d s\n", "Remaining CPU time", day, hour, min, sec);

    if ((rest < 86400)&&(rest >= 3600))
      fprintf(f, "%23s : %2d h %2d min %2d s\n", "Remaining CPU time", hour, min, sec);

    if ((rest < 3600)&&(rest >= 60))
      fprintf(f, "%23s : %2d min %2d s\n", "Remaining CPU time", min, sec);

    if (rest < 60)
      fprintf(f, "%23s : %2d s\n", "Remaining CPU time", sec);
  }
      
  if (Verbose)
  {
    printf("\033[1A\033[1A\033[1A");	     
    printf("Iteration No %i\n", IterationNo);
    printf("Time step = %15.8e\n",TimeStep);
    printf("Computation at %5.1f %1s \n",(TotalTime-RemainingTime)*100./TotalTime,"%");    
  }
  
  fprintf(f,"\n");
  
  fclose(f);
  
}

/*
_______________________________________________________________________________

Indicate the time to write
_______________________________________________________________________________

*/

bool timeToWrite()
{
  return ((TotalTime - RemainingTime) > (TotalTime*WriteNo/WriteNb));
}
