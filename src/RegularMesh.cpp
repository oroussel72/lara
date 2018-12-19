#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
using namespace std;

#include "Global.h"
#include "Vector.h"
#include "Matrix.h"
#include "Timer.h"
#include "Parameters.h"
#include "Cell.h"
#include "RegularMesh.h"

/*
_______________________________________________________________________________

Constructor
_______________________________________________________________________________

*/

RegularMesh::RegularMesh()
{
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  int i,j,k;
  
  real dx = (XMax - XMin)/NX;
  real dy = (YMax - YMin)/NY;
  real dz = (ZMax - ZMin)/NZ;
    
  // Init array of cells
  
  switch (Dimension)
  {
    case 1:
      _cell = new Cell[NX+6];
      break;
    case 2:
      _cell = new Cell[(NX+6)*(NY+6)];
      break;
    case 3:
      _cell = new Cell[(NX+6)*(NY+6)*(NZ+6)];
      break;
  };
  
  // Set mesh
  
#pragma omp parallel for private(i,j,k) shared(dx,dy,dz,ei,ej,ek)
  for (i=0; i<=ei*(NX+5); i++)
  for (j=0; j<=ej*(NY+5); j++)
  for (k=0; k<=ek*(NZ+5); k++)
  {
    cell(i,j,k)->setCenter(1,XMin + (i-3+0.5)*dx);
    if (Dimension > 1) cell(i,j,k)->setCenter(2,YMin + (j-3+0.5)*dy);
    if (Dimension > 2) cell(i,j,k)->setCenter(3,ZMin + (k-3+0.5)*dz);
    
    cell(i,j,k)->setSize(1,dx);
    if (Dimension > 1) cell(i,j,k)->setSize(2,dy);
    if (Dimension > 2) cell(i,j,k)->setSize(3,dz);
  }
  
  // Init space step
  
  SpaceStep = dx;
  if (Dimension > 1) SpaceStep = min(SpaceStep,dy);
  if (Dimension > 2) SpaceStep = min(SpaceStep,dz);
  
}

/*
_______________________________________________________________________________

Distructor
_______________________________________________________________________________

*/

RegularMesh::~RegularMesh()
{
  delete[] _cell;
}

/*
_______________________________________________________________________________

Compute divergence
_______________________________________________________________________________

*/

void RegularMesh::computeDivergence()
{
  int i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  Vector Fx1,Fx2,Fy1,Fy2,Fz1,Fz2;
  
  // Loop on real cells
  
#pragma omp parallel for private(i,j,k,Fx1,Fx2,Fy1,Fy2,Fz1,Fz2) shared(ei,ej,ek)
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {
      // Compute fluxes in the x-direction
    
      Fx2 = flux(cell(i-2,j,k),cell(i-1,j,k),cell(i,j,k),cell(i+1,j,k),cell(i+2,j,k),cell(i+3,j,k),1);
      Fx1 = flux(cell(i-3,j,k),cell(i-2,j,k),cell(i-1,j,k),cell(i,j,k),cell(i+1,j,k),cell(i+2,j,k),1);
    
      cell(i,j,k)->setDivergence((Fx1-Fx2)/cell(i,j,k)->dx());
    
      if (Dimension > 1)
      {
        // Compute fluxes in the y-direction
      
        Fy2 = flux(cell(i,j-2,k),cell(i,j-1,k),cell(i,j,k),cell(i,j+1,k),cell(i,j+2,k),cell(i,j+3,k),2);
        Fy1 = flux(cell(i,j-3,k),cell(i,j-2,k),cell(i,j-1,k),cell(i,j,k),cell(i,j+1,k),cell(i,j+2,k),2);
    
        cell(i,j,k)->addDivergence((Fy1-Fy2)/cell(i,j,k)->dy());   
      }

      if (Dimension > 2)
      {
        // Compute fluxes in the z-direction
      
        Fz2 = flux(cell(i,j,k-2),cell(i,j,k-1),cell(i,j,k),cell(i,j,k+1),cell(i,j,k+2),cell(i,j,k+3),3);
        Fz1 = flux(cell(i,j,k-3),cell(i,j,k-2),cell(i,j,k-1),cell(i,j,k),cell(i,j,k+1),cell(i,j,k+2),3);
    
        cell(i,j,k)->addDivergence((Fz1-Fz2)/cell(i,j,k)->dz()); 
      }  
    
  }
}

/*
_______________________________________________________________________________

Check the stability of the computation
_______________________________________________________________________________

*/

void RegularMesh::checkStability()
{
  int i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  FILE* f;
  
  // Loop on real cells
  
#pragma omp parallel for private(i,j,k) shared(ei,ej,ek)
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {
    if (cell(i,j,k)->isNaN())
    {
      printf("lara: Numerical instability detected at iteration %i, position (%15.8f,%15.8f,%15.8f)\n", IterationNo, cell(i,j,k)->x(), (Dimension > 1) ? cell(i,j,k)->y():0.,
	     (Dimension > 2) ? cell(i,j,k)->z():0.);
      f=fopen("performance.out","a");
      fprintf(f,"Numerical instability detected at iteration %i, position (%15.8f,%15.8f,%15.8f)\n\n", IterationNo, cell(i,j,k)->x(), (Dimension > 1) ? cell(i,j,k)->y():0.,
	     (Dimension > 2) ? cell(i,j,k)->z():0.);
      fclose(f);
      exit(-1);
    }
  }
}

/*
_______________________________________________________________________________

Copy cell-average values from sourceStageNo to targetStageNo
_______________________________________________________________________________

*/

void RegularMesh::copyAverage(const int sourceStorageNo, const int targetStorageNo)
{
  int i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  // Loop on real cells
  
#pragma omp parallel for private(i,j,k) shared(ei,ej,ek)
  for (i=ei; i <= ei*(NX + 5); i++)
  for (j=ej; j <= ej*(NY + 5); j++)
  for (k=ek; k <= ek*(NZ + 5); k++)
    cell(i,j,k)->copyStorage(sourceStorageNo, targetStorageNo);
} 

/*
_______________________________________________________________________________

Compute one step of a multi-stage time integration method
_______________________________________________________________________________

*/

void RegularMesh::computeStep()
{
  int i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  real X,Y,Z,t;
  
  real a = 1./3.;
   
  Vector U0, U, D;
  real dt=TimeStep;
  
  X=Y=Z=t=0.;
  
  // Loop on real cells
  
#pragma omp parallel for private(i,j,k,U,U0,D) shared(dt,ei,ej,ek,a)
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {    
      switch (StepNo)
      {
        case 1:
        default:	
	  U = cell(i,j,k)->average();
	  D = cell(i,j,k)->divergence();
	  cell(i,j,k)->setAverage(U + dt*D);
	  break;

        case 2:
	  switch(StepNb)
	  {
	    case 1:
	    default:
	      break;
      
	    case 2:
	      U0 = cell(i,j,k)->average(2);
	      U  = cell(i,j,k)->average();
	      D  = cell(i,j,k)->divergence();
	      cell(i,j,k)->setAverage(0.5*(U0 + U + dt*D) );
	      break;
	    
	    case 3:
	      U0 = cell(i,j,k)->average(2);
	      U  = cell(i,j,k)->average();
	      D  = cell(i,j,k)->divergence();
	      cell(i,j,k)->setAverage(0.25*(3.*U0 + U + dt*D) );
	      break;
	  };
	  break;
	
        case 3:
	  switch(StepNb)
	  {
	    case 1:
	    default:
	      break;
	      
	    case 2:
	      break;
	      
	    case 3:
	      U0 = cell(i,j,k)->average(2);
	      U  = cell(i,j,k)->average();
	      D  = cell(i,j,k)->divergence();
	      cell(i,j,k)->setAverage(a*(U0 + 2.*(U + dt*D)) );
	      break;    
  	  };
	  break;
      }; 
      
      if (UsePenalization)
      {
	X = cell(i,j,k)->x();
	Y = cell(i,j,k)->y();
	Z = cell(i,j,k)->z();
	t = TotalTime-RemainingTime;
	cell(i,j,k)->setAverage(QuantityNb,1,(isImmersed(X,Y,Z,t)? 0.:1.));
      }
  }
}

/*
_______________________________________________________________________________

Compute integral values
_______________________________________________________________________________

*/

void RegularMesh::computeIntegralValues()
{
  int i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  real dx, dy, dz;
  real x, y, z;
  int axisNo,speciesNo;
  real mu;
  
  MaxEigenvalue = 0.;
  KineticEnergy = 0.;
  MaxDiffusivity = (Equation==NavierStokes)? 0.: Diffusivity;
  MaxAverage.reset();
  
  if (IterationNo > 0)
    AverageTimeStep = ((IterationNo-1)*AverageTimeStep + TimeStep)/IterationNo;
  else
    AverageTimeStep = TimeStep;
  
  // Loop on real cells
  
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {  
    MaxEigenvalue = max(MaxEigenvalue,cell(i,j,k)->maxEigenvalue());
    
    // --- Compute maximal diffusivity ---

    if (Equation == NavierStokes)
    {
      mu = cell(i,j,k)->sutherland(cell(i,j,k)->T());
      MaxDiffusivity = MAX(MaxDiffusivity, MAX(Pr,1.)*mu/cell(i,j,k)->rho());    
    }  
      
    dx = cell(i,j,k)->dx();
    dy = (Dimension > 1)? cell(i,j,k)->dy():1.;
    dz = (Dimension > 2)? cell(i,j,k)->dz():1.;
    x = cell(i,j,k)->x();
    y = (Dimension > 1) ? cell(i,j,k)->y():0.;
    z = (Dimension > 2) ? cell(i,j,k)->z():0.;
    
    switch(Equation)
    { 
      case ConvectionDiffusion:
	KineticEnergy += ABS(cell(i,j,k)->density()-initCondition(x,y,z).at(1))*dx*dy*dz;
	break;
	
      case Burgers:
	KineticEnergy += 0.5*SQR(cell(i,j,k)->velocityNorm())*dx*dy*dz;
	break;
	
      case NavierStokes:
	KineticEnergy += 0.5*cell(i,j,k)->density()*SQR(cell(i,j,k)->velocityNorm())*dx*dy*dz;
	MaxAverage.set(1,MAX(MaxAverage.at(1),ABS(cell(i,j,k)->density())));
	for (axisNo=1; axisNo<=Dimension; axisNo++)
	  MaxAverage.set(axisNo+1,MAX(MaxAverage.at(axisNo+1),cell(i,j,k)->momentumNorm()));
	MaxAverage.set(Dimension+2,MAX(MaxAverage.at(Dimension+2),ABS(cell(i,j,k)->energy())));
	if (SpeciesNb > 1)
	{
	  for (speciesNo = 1; speciesNo < SpeciesNb; speciesNo++)
	    MaxAverage.set(Dimension+2+speciesNo,MAX(MaxAverage.at(Dimension+2+speciesNo),cell(i,j,k)->average(Dimension+2+speciesNo,1)));
	}
	break;
    };
  }
  
  // Compute maximum of the density averaged in time (for detonation only)
  
  if (Detonation)
  {
    if (IterationNo > 0)
      TimeAverageMaxDensity = ((IterationNo-1.)*TimeAverageMaxDensity+MaxAverage.at(1))/(1.*IterationNo);
    else
      TimeAverageMaxDensity = 0.;
  }
}

/*
_______________________________________________________________________________

Compute time step in function of CFL
_______________________________________________________________________________

*/

void RegularMesh::computeTimeStep()
{
  if (CFL != 0.)
  {

    if (IsInviscid)
      TimeStep = CFL*SpaceStep/MaxEigenvalue;
    else
      TimeStep = CFL*MIN(SpaceStep/MaxEigenvalue, SQR(SpaceStep)/(2.*MaxDiffusivity));
  
    if (TimeStep > RemainingTime) TimeStep = RemainingTime;
  }
    
}

/*
_______________________________________________________________________________

Compute gradient
_______________________________________________________________________________

*/

void RegularMesh::computeGradient()
{
  int i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  // Loop on real cells
  
#pragma omp parallel for private(i,j,k) shared(ei,ej,ek)
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {
    cell(i,j,k)->computeGradient(cell(i-2,j,k),cell(i-1,j,k),cell(i+1,j,k),cell(i+2,j,k),1);
    if (Dimension > 1) cell(i,j,k)->computeGradient(cell(i,j-2,k),cell(i,j-1,k),cell(i,j+1,k),cell(i,j+2,k),2);
    if (Dimension > 2) cell(i,j,k)->computeGradient(cell(i,j,k-2),cell(i,j,k-1),cell(i,j,k+1),cell(i,j,k+2),3);
  }
}

/*
_______________________________________________________________________________

Compute source
_______________________________________________________________________________

*/

void RegularMesh::computeSource()
{
  int i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  // Loop on real cells
  
#pragma omp parallel for private(i,j,k) shared(ei,ej,ek)
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
    cell(i,j,k)->computeSource();
}
/*
_______________________________________________________________________________

Set boundary conditions
_______________________________________________________________________________

*/

void RegularMesh::setBoundaryConditions()
{
  int i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  int imin = ei*3;
  int jmin = ej*3;
  int kmin = ek*3;
  int imax = ei*(NX+2);
  int jmax = ej*(NY+2);
  int kmax = ek*(NZ+2);
    
    // --- XMin ---------------------------------------------------------------

#pragma omp parallel for private(j,k) shared(jmin,jmax,kmin,kmax)
    for (j=jmin; j <= jmax; j++)
    for (k=kmin; k <= kmax; k++)
    {
	cell(0,j,k)->linkBoundaryCell(BCXMin,cell(5,j,k),cell(1,j,k),cell(NX  ,j,k),1);
	cell(1,j,k)->linkBoundaryCell(BCXMin,cell(4,j,k),cell(2,j,k),cell(NX+1,j,k),1);
	cell(2,j,k)->linkBoundaryCell(BCXMin,cell(3,j,k),cell(3,j,k),cell(NX+2,j,k),1);
    }
    
    // --- XMax ---------------------------------------------------------------
        
#pragma omp parallel for private(j,k) shared(jmin,jmax,kmin,kmax)
    for (j=jmin; j <= jmax; j++)
    for (k=kmin; k <= kmax; k++)
    {
	  cell(NX+3,j,k)->linkBoundaryCell(BCXMax,cell(NX+2,j,k),cell(NX+2,j,k),cell(3,j,k),1);
	  cell(NX+4,j,k)->linkBoundaryCell(BCXMax,cell(NX+1,j,k),cell(NX+3,j,k),cell(4,j,k),1);
	  cell(NX+5,j,k)->linkBoundaryCell(BCXMax,cell(NX  ,j,k),cell(NX+4,j,k),cell(5,j,k),1);
    }
 
  
    if (Dimension > 1)
    {
      // --- YMin -------------------------------------------------------------
      
#pragma omp parallel for private(i,k) shared(imin,imax,kmin,kmax)
	for (i=imin; i <= imax; i++)
	for (k=kmin; k <= kmax; k++)
	{
	    cell(i,0,k)->linkBoundaryCell(BCYMin,cell(i,5,k),cell(i,1,k),cell(i,NY  ,k),2);
	    cell(i,1,k)->linkBoundaryCell(BCYMin,cell(i,4,k),cell(i,2,k),cell(i,NY+1,k),2);
	    cell(i,2,k)->linkBoundaryCell(BCYMin,cell(i,3,k),cell(i,3,k),cell(i,NY+2,k),2);
	}
    
      
      // --- YMax -------------------------------------------------------------
      
#pragma omp parallel for private(i,k) shared(imin,imax,kmin,kmax)
	for (i=imin; i <= imax; i++)
	for (k=kmin; k <= kmax; k++)
	{
	    cell(i,NY+3,k)->linkBoundaryCell(BCYMax,cell(i,NY+2,k),cell(i,NY+2,k),cell(i,3,k),2);
	    cell(i,NY+4,k)->linkBoundaryCell(BCYMax,cell(i,NY+1,k),cell(i,NY+3,k),cell(i,4,k),2);
	    cell(i,NY+5,k)->linkBoundaryCell(BCYMax,cell(i,NY  ,k),cell(i,NY+4,k),cell(i,5,k),2);
	}
      
    }
   
    if (Dimension > 2)
    {
      
      // --- ZMin -------------------------------------------------------------
      
#pragma omp parallel for private(i,j) shared(imin,imax,jmin,jmax)
	for (i=imin; i <= imax; i++)
	for (j=jmin; j <= jmax; j++)
	{
	    cell(i,j,0)->linkBoundaryCell(BCZMin,cell(i,j,5),cell(i,j,1),cell(i,j,NZ  ),3);
	    cell(i,j,1)->linkBoundaryCell(BCZMin,cell(i,j,4),cell(i,j,2),cell(i,j,NZ+1),3);
	    cell(i,j,2)->linkBoundaryCell(BCZMin,cell(i,j,3),cell(i,j,3),cell(i,j,NZ+2),3);
	}
      
      // --- ZMax -------------------------------------------------------------
      
#pragma omp parallel for private(i,j) shared(imin,imax,jmin,jmax)
	for (i=imin; i <= imax; i++)
	for (j=jmin; j <= jmax; j++)
	{
	    cell(i,j,NZ+3)->linkBoundaryCell(BCZMax,cell(i,j,NZ+2),cell(i,j,NZ+2),cell(i,j,3),3);
	    cell(i,j,NZ+4)->linkBoundaryCell(BCZMax,cell(i,j,NZ+1),cell(i,j,NZ+3),cell(i,j,4),3);
	    cell(i,j,NZ+5)->linkBoundaryCell(BCZMax,cell(i,j,NZ  ),cell(i,j,NZ+4),cell(i,j,5),3);
	}
      
    }    
}

/*
_______________________________________________________________________________

Compute complete iteration
_______________________________________________________________________________

*/

void RegularMesh::iterate()
{
  // Init iteration
  
  computeTimeStep();

  // Set boundary conditions 
  
  setBoundaryConditions();
    
  // For detonation computations, compute source term on Dt/2
  
  if (Detonation) computeSource();
  
  if (UsePenalization) 
  { 
    render();
    penalize();
  }
  
  // Multi-stage method
  
  copyAverage(1,2);
    
  for (StepNo = 1; StepNo <= StepNb; StepNo++)
  {
    if (Equation == NavierStokes && !IsInviscid) 
      computeGradient();
    
    computeDivergence();
    computeStep();
    checkStability();
  }

  // For detonation computations, compute source term on Dt/2
  
  if (Detonation) computeSource();

  if (UsePenalization) penalize();
  
  // Store old average values
  
  if (IterationNo > 1) copyAverage(3,4);
  copyAverage(1,3);
  
  // Compute integral values
  
  computeIntegralValues();

  // Write cell-average values into file if required
  
  if (timeToWrite())
  {
     writeAverage();
     backup();
  }
  
  // If it is time, show the evolution of the computation
     
  if (timeToShow())  
  {
    showPerformance();
    writeIntegral();
  }  
  
  // Compute remaining time
  
  RemainingTime -= TimeStep;
}

/*
_______________________________________________________________________________

Write cell-average values
_______________________________________________________________________________

*/

void RegularMesh::writeAverage()
{
  int i,j,k,l;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  int imin = ei*3;
  int imax = ei*(NX+2);
  int jmin = ej*3;
  int jmax = ej*(NY+2);
  int kmin = ek*3;
  int kmax = ek*(NZ+2);
  
  real dx = (XMax-XMin)/NX;
  real dy = (YMax-YMin)/NY;
  real dz = (ZMax-ZMin)/NZ;
    
  char fileName[256];
  FILE *f;
  
  // Write data
  
  if (WriteFormat == VTK)
    sprintf(fileName,"Average%05i.vtk",WriteNo);
  else
    sprintf(fileName,"Average%05i.gnu",WriteNo);
  
  // Open file
  
  f = fopen(fileName,"w");
    
  switch (WriteFormat)
  {  
    case VTK:
      
      // --- HEADER ---
      
      fprintf(f,"# vtk DataFile Version 3.0\n");
      fprintf(f,"Generated by Lara\n");
      
      if (WriteType == ASCII)
	fprintf(f,"ASCII\n");
      else
	fprintf(f,"BINARY\n");
	      
      fprintf(f,"DATASET STRUCTURED_POINTS\n");
      fprintf(f,"DIMENSIONS %i %i %i\n", NX+1, (Dimension < 2)? 2:NY+1, (Dimension < 3)? 2:NZ+1);
      fprintf(f,"ORIGIN %f %f %f\n", XMin*ei, YMin*ej, ZMin*ek);
      fprintf(f,"SPACING %f %f %f\n", dx, (Dimension < 2)? dx:dy, (Dimension < 3) ? dx:dz);
      fprintf(f,"CELL_DATA %i\n", NX*NY*NZ);

      // --- DATA ---
      
      if (Equation != Burgers)
      {
	if (Equation == ConvectionDiffusion)
	  fprintf(f,"SCALARS u float\n");
	else
	  fprintf(f,"SCALARS rho float\n");
	
	fprintf(f,"LOOKUP_TABLE default\n");    
	
	for (k=kmin; k<=kmax; k++)
	for (j=jmin; j<=jmax; j++)
	for (i=imin; i<=imax; i++)
	  writeLn(f,cell(i,j,k)->rho());
      }    
      
      if (Equation == NavierStokes)
      {
	if (WriteType==Binary) fprintf(f,"\n");
	fprintf(f,"SCALARS p float\n");
	fprintf(f,"LOOKUP_TABLE default\n");
      
	for (k=kmin; k<=kmax; k++)
	for (j=jmin; j<=jmax; j++)
	for (i=imin; i<=imax; i++)
	  writeLn(f,cell(i,j,k)->p());
      
	if (WriteType==Binary) fprintf(f,"\n");
	fprintf(f,"SCALARS T float\n");
	fprintf(f,"LOOKUP_TABLE default\n");
      
	for (k=kmin; k<=kmax; k++)
	for (j=jmin; j<=jmax; j++)
	for (i=imin; i<=imax; i++)
	  writeLn(f,cell(i,j,k)->T());
      }    
      
      if (Equation != ConvectionDiffusion)
      {
	if (Dimension == 1)
	{
	  if (WriteType==Binary) fprintf(f,"\n");
	  fprintf(f,"SCALARS u float\n");
	  fprintf(f,"LOOKUP_TABLE default\n");
      
	  for (k=kmin; k<=kmax; k++)
	  for (j=jmin; j<=jmax; j++)
	  for (i=imin; i<=imax; i++)
	    writeLn(f,cell(i,j,k)->v(1));
	}
	else
	{
	  if (WriteType==Binary) fprintf(f,"\n");
	  fprintf(f,"VECTORS v float\n");
	
	  for (k=kmin; k<=kmax; k++)
	  for (j=jmin; j<=jmax; j++)
	  for (i=imin; i<=imax; i++)
	  {
	    write(f,cell(i,j,k)->v(1));
	    write(f,cell(i,j,k)->v(2));
	    write(f,(Dimension > 2) ? cell(i,j,k)->v(3):0.);
	    if (WriteType==ASCII) fprintf(f,"\n"); 
	  }
	}
      }
      
      if (Equation == NavierStokes && SpeciesNb > 1)
      {
	for (l=1; l<=SpeciesNb-1; l++)
	{
	  if (WriteType==Binary) fprintf(f,"\n");
	  
	  if (UsePenalization && l==SpeciesNb-1)
	    fprintf(f,"SCALARS Mask float\n");
	  else
	    fprintf(f,"SCALARS Y%i float\n",l);
	  
	  fprintf(f,"LOOKUP_TABLE default\n");
	  
	  for (k=kmin; k<=kmax; k++)
	  for (j=jmin; j<=jmax; j++)
	  for (i=imin; i<=imax; i++)
	  {
	    if (UsePenalization && l==SpeciesNb-1)
	      writeLn(f,cell(i,j,k)->average(QuantityNb,1));
	    else
	      writeLn(f,cell(i,j,k)->Y(l));
	  }
	} 
      }
      
      fprintf(f,"\n");
    break;
    
  case Gnuplot:
   
    // --- HEADER ---
    
    fprintf(f, "# %13s ", "x");
    if (Dimension > 1) fprintf(f, "%15s ", "y");
    if (Dimension > 2) fprintf(f, "%15s ", "z");
    if (Equation == ConvectionDiffusion) fprintf(f, "%15s ", "u");
    if (Equation == NavierStokes) fprintf(f, "%15s %15s %15s ", "rho", "p", "T");

    if (Equation != ConvectionDiffusion)
    {
      fprintf(f, "%15s ", "u");
      if (Dimension > 1) fprintf(f, "%15s ", "v");
      if (Dimension > 2) fprintf(f, "%15s ", "w");
    }
    
    if (Equation == NavierStokes && SpeciesNb > 1)
    {
      for (l=1; l<=SpeciesNb-1; l++)
	fprintf(f,"             Y%i ",l); 
    }
    
    fprintf(f,"\n");

    // --- DATA ---
    
    for (k=kmin; k<=kmax; k++)
    {
      for (j=jmin; j<=jmax; j++)
      {
	for (i=imin; i<=imax; i++)
	{
	  for (l=1; l<=Dimension;l++)
	    write(f,cell(i,j,k)->center(l));
	      
	  if (Equation != Burgers)    
	    write(f,cell(i,j,k)->rho());
	  
	  if (Equation == NavierStokes)
	  {
	    write(f,cell(i,j,k)->p());
	    write(f,cell(i,j,k)->T());
	  }
	  
	  if (Equation != ConvectionDiffusion)
	  {
	    for (l=1; l<=Dimension;l++)
	      write(f,cell(i,j,k)->v(l));
	  }
	  
	  if (Equation == NavierStokes && SpeciesNb > 1)
	  {
	    for (l=1; l<=SpeciesNb-1;l++)
	      write(f,cell(i,j,k)->Y(l));
	  }
	      
	  fprintf(f,"\n");
	}
	fprintf(f,"\n");	
      }
      fprintf(f,"\n");
    }
    break;  
  };
	  
  fclose(f);
  
  WriteNo++;
}

/*
_______________________________________________________________________________

Init mesh
_______________________________________________________________________________

*/

void RegularMesh::init()
{
  cout << "lara: execution\n";
  
  initAverage();
  if (Restart) restore();
  copyAverage(1,3);
  copyAverage(1,4);
  computeIntegralValues();

  if (!Restart)
  {
    writeAverage();
    writeIntegral();
  }
  
  if (Verbose)
    cout << endl << endl << endl;
  else
    cout << "lara: task launch in background\n";
}

/*
_______________________________________________________________________________

End mesh
_______________________________________________________________________________

*/

void RegularMesh::end()
{
  writeAverage();
  computeIntegralValues();
  writeIntegral();
  showPerformance();
}

/*
_______________________________________________________________________________

Compute initial condition
_______________________________________________________________________________

*/

void RegularMesh::initAverage()
{
  int i,j,k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  real x,y,z;
  
#pragma omp parallel for private(i,j,k,x,y,z) shared(ei,ej,ek)
  for (i=0; i<=ei*(NX+5); i++)
  for (j=0; j<=ej*(NY+5); j++)
  for (k=0; k<=ek*(NZ+5); k++)
  {
    x = cell(i,j,k)->x();
    y = (Dimension > 1) ? cell(i,j,k)->y() : 0.;
    z = (Dimension > 2) ? cell(i,j,k)->z() : 0.;
    cell(i,j,k)->setAverage(1,initCondition(x,y,z));
  }
}

/*
_______________________________________________________________________________

Write integral values
_______________________________________________________________________________

*/

void RegularMesh::writeIntegral()
{
  FILE *f;
  
  real delta = sqrt(2*(Gamma-1.)/(Gamma+1.));
  real rhoVN = Gamma/(1-delta);
  real err = ABS(rhoVN-TimeAverageMaxDensity)/ABS(rhoVN);
  
  
  if (IterationNo == 0 && !Restart)
  {  
    f = fopen("Integral.gnu","w");
    fprintf(f,"# %13s %15s %15s %15s\n", "Time", "Time step", "Kin. energy", "Max density");
  }
  else
    f = fopen("Integral.gnu","a");
  
  fprintf(f,"%15.8e %15.8e %15.8e %15.8e\n", TotalTime-RemainingTime, TimeStep, KineticEnergy, ((Detonation)? err:MaxAverage.at(1)));
  
  fclose(f);
}

/*
_______________________________________________________________________________

Backup data
_______________________________________________________________________________

*/

void RegularMesh::backup()
{
  FILE *f;
  int i,j,k,n;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  f = fopen("lara.bak","w");
  
  fprintf(f,"%15.8e\n", TimeStep);
  fprintf(f,"%15.8e\n", TotalTime-RemainingTime);
  fprintf(f,"%15i\n", IterationNo);
  fprintf(f,"%15i\n", WriteNo);
  
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {
    for (n=1;n<=QuantityNb;n++)
      fprintf(f,"%15.8e  ",cell(i,j,k)->average(n,1));
    
    fprintf(f,"\n");
  }
  fclose(f);
}

/*
_______________________________________________________________________________

Restore data
_______________________________________________________________________________

*/

void RegularMesh::restore()
{
  FILE *f;
  int i,j,k,n;
  int ei = 1;
  int rep = 0;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  float elapsedTime, timeStep;
  float x;
  
  f = fopen("lara.bak","r");
  
  if (!f) return;
  
  rep = fscanf(f,"%e", &timeStep);
  rep = fscanf(f,"%e", &elapsedTime);
  rep = fscanf(f,"%i", &IterationNo);
  rep = fscanf(f,"%i", &WriteNo);
 
  TimeStep = timeStep;
  RemainingTime = TotalTime - elapsedTime;
  StartTime = elapsedTime;
  
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {
    for (n=1;n<=QuantityNb;n++)
    {
      rep = fscanf(f,"%e  ",&x);
      cell(i,j,k)->setAverage(n,1,x);
    }
  }
  
  fclose(f);
}

/*
_______________________________________________________________________________

Penalize 
_______________________________________________________________________________

*/

void RegularMesh::penalize()
{
  int i,j,k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;

// --- penalize momentum n->n+1/2 ---
  
#pragma omp parallel for private(i,j,k) shared(ei,ej,ek)
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {
    if (cell(i,j,k)->isInSolid()) 
      cell(i,j,k)->penalizeMomentum();
  }
  
// --- penalize density n->n+1 ---
  
 #pragma omp parallel for private(i,j,k) shared(ei,ej,ek)
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {
    if (cell(i,j,k)->isInSolid()) 
    {
	cell(i,j,k)->setDensity(cell(i,j,k)->density()+TimeStep/(4.*cell(i,j,k)->size(1))*(cell(i-1,j,k)->momentum(1)-cell(i+1,j,k)->momentum(1)));

	if (Dimension > 1)
	  cell(i,j,k)->setDensity(cell(i,j,k)->density()+TimeStep/(4.*cell(i,j,k)->size(2))*(cell(i,j-1,k)->momentum(2)-cell(i,j+1,k)->momentum(2)));	  
	
	if (Dimension > 2)
	  cell(i,j,k)->setDensity(cell(i,j,k)->density()+TimeStep/(4.*cell(i,j,k)->size(3))*(cell(i,j,k-1)->momentum(3)-cell(i,j,k+1)->momentum(3)));	  
	
/*
	cell(i,j,k)->penalizeDensity(cell(i-3,j,k), cell(i-2,j,k),cell(i-1,j,k), cell(i+1,j,k),cell(i+2,j,k), cell(i+3,j,k),1);

	if (Dimension > 1) 
	  cell(i,j,k)->penalizeDensity(cell(i,j-3,k), cell(i,j-2,k),cell(i,j-1,k), cell(i,j+1,k),cell(i,j+2,k), cell(i,j+3,k),2);
	
	if (Dimension > 2) 
	  cell(i,j,k)->penalizeDensity(cell(i,j,k-3), cell(i,j,k-2),cell(i,j,k-1), cell(i,j,k+1),cell(i,j,k+2), cell(i,j,k+3),3);
*/	
    }
  }   
    
// --- penalize momentum n+1/2->n+1 ---
  
#pragma omp parallel for private(i,j,k) shared(ei,ej,ek)
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {
    if (cell(i,j,k)->isInSolid()) 
	cell(i,j,k)->penalizeMomentum();
  }
  
// --- penalize momentum n->n+1 ---
  
#pragma omp parallel for private(i,j,k) shared(ei,ej,ek)
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {
    if (cell(i,j,k)->isInSolid())   
	cell(i,j,k)->penalizeEnergy();
  }
}
/*
_______________________________________________________________________________

Interpolate in cell l,i,j,k when the cell becomes fluid
_______________________________________________________________________________

*/

void RegularMesh::interpolateFluid(const int i, const int j, const int k)
{
  int n;
  Vector result(QuantityNb);
  
  n = 0;
  
  if (!cell(i-1,j,k)->isInSolid())
  {
    n++;
    result += cell(i-1,j,k)->average();
  }
  
  if (!cell(i+1,j,k)->isInSolid())
  {
    n++;
    result += cell(i+1,j,k)->average();
  }
  
  if (Dimension > 1 && !cell(i,j-1,k)->isInSolid())
  {
    n++;
    result += cell(i,j-1,k)->average();
  }
     
  if (Dimension > 1 && !cell(i,j+1,k)->isInSolid())
  {
    n++;
    result += cell(i,j+1,k)->average();
  }
     
   if (Dimension > 2 && !cell(i,j,k-1)->isInSolid())
  {
    n++;
    result += cell(i,j,k-1)->average();
  }
     
  if (Dimension > 1 && !cell(i,j,k+1)->isInSolid())
  {
    n++;
    result += cell(i,j,k+1)->average();
  }
     
  result *= (1./n);
  
  cell(i,j,k)->setAverage(result);
  printf("interpolateFluid\n");
}


/*
_______________________________________________________________________________

Interpolate in cell l,i,j,k when the cell becomes solid
_______________________________________________________________________________

*/

void RegularMesh::interpolateSolid(const int i, const int j, const int k)
{
  int n;
  Vector result(QuantityNb);
  
  n = 0;
  
  if (cell(i-1,j,k)->isInSolid())
  {
    n++;
    result += cell(i-1,j,k)->average();
  }
  
  if (cell(i+1,j,k)->isInSolid())
  {
    n++;
    result += cell(i+1,j,k)->average();
  }
  
  if (Dimension > 1 && cell(i,j-1,k)->isInSolid())
  {
    n++;
    result += cell(i,j-1,k)->average();
  }
     
  if (Dimension > 1 && cell(i,j+1,k)->isInSolid())
  {
    n++;
    result += cell(i,j+1,k)->average();
  }
     
   if (Dimension > 2 && cell(i,j,k-1)->isInSolid())
  {
    n++;
    result += cell(i,j,k-1)->average();
  }
     
  if (Dimension > 1 && cell(i,j,k+1)->isInSolid())
  {
    n++;
    result += cell(i,j,k+1)->average();
  }
     
  result *= (1./n);
  
  cell(i,j,k)->setAverage(result);
}

/*
_______________________________________________________________________________

Render date
_______________________________________________________________________________

*/
void RegularMesh::render()
{
  int i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  if (PenalVelocity == 0.) return;
  
#pragma omp parallel for private(i,j,k) shared(ei,ej,ek)
  for (i=ei*3; i <= ei*(NX + 2); i++)
  for (j=ej*3; j <= ej*(NY + 2); j++)
  for (k=ek*3; k <= ek*(NZ + 2); k++)
  {    
    if (cell(i,j,k)->status()==Leaf)
    {
      if(cell(i,j,k)->wasInSolid() && !cell(i,j,k)->isInSolid())
	interpolateFluid(i,j,k);
      
      if(!cell(i,j,k)->wasInSolid() && cell(i,j,k)->isInSolid())
	interpolateSolid(i,j,k);
    }
  }
}

