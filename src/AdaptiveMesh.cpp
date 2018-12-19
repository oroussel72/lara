#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
using namespace std;

#include "Global.h"
#include "Vector.h"
#include "Matrix.h"
#include "Timer.h"
#include "Parameters.h"
#include "Cell.h"
#include "AdaptiveMesh.h"

/*
_______________________________________________________________________________

Constructor
_______________________________________________________________________________

*/

AdaptiveMesh::AdaptiveMesh()
{
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  int l,i,j,k;
  
  real dx,dy,dz;
  
  dx = dy = dz = 0.;
    
  // Init S
  
  _cellNb = new int[ScaleNb+2];
  
  if (Dimension == 1)  
    for (l=0; l<=ScaleNb+1; l++) _cellNb[l] = (1<<l) + 6*l - 1;
  else if (Dimension == 2)
    for (l=0; l<=ScaleNb+1; l++) _cellNb[l] = ((1<<(2*l))-1)/3 + 12*(1<<l) + 36*l - 12;
  else
    for (l=0; l<=ScaleNb+1; l++) _cellNb[l] = ((1<<(3*l))-1)/7 + 6*(1<<(2*l)) + 108*(1<<l) + 216*l - 114;
  
  // Init array of cells
  
  _cell = new Cell[_cellNb[ScaleNb+1]];
  
  // Set mesh
  
  for (l=0; l<=ScaleNb; l++)
  {  
    dx = (XMax-XMin)/N(l);
    dy = (YMax-YMin)/N(l);
    dz = (ZMax-ZMin)/N(l);

#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek,dx,dy,dz)
    for (i=-3*ei; i<=ei*(N(l)+2); i++)
    for (j=-3*ej; j<=ej*(N(l)+2); j++)
    for (k=-3*ek; k<=ek*(N(l)+2); k++)
    {
      cell(l,i,j,k)->setCenter(1,XMin + (i+0.5)*dx);
      if (Dimension > 1) cell(l,i,j,k)->setCenter(2,YMin + (j+0.5)*dy);
      if (Dimension > 2) cell(l,i,j,k)->setCenter(3,ZMin + (k+0.5)*dz);
    
      cell(l,i,j,k)->setSize(1,dx);
      if (Dimension > 1) cell(l,i,j,k)->setSize(2,dy);
      if (Dimension > 2) cell(l,i,j,k)->setSize(3,dz);
    }
  }
  
  // Init space step
  
  SpaceStep = dx;
  if (Dimension > 1) SpaceStep = min(SpaceStep,dy);
  if (Dimension > 2) SpaceStep = min(SpaceStep,dz);
  
}

/*
_______________________________________________________________________________

Constructor
_______________________________________________________________________________

*/

AdaptiveMesh::~AdaptiveMesh()
{
  delete[] _cell;
  delete[] _cellNb;
}

/*
_______________________________________________________________________________

Split cell
_______________________________________________________________________________

*/

void AdaptiveMesh::split(const int l, const int i, const int j, const int k)
{
  int ni,nj,nk;
  
  int ei = 1;
  int ej = (Dimension > 1)? 1:0;
  int ek = (Dimension > 2)? 1:0;
  
  cell(l,i,j,k)->setStatus(Node);
  
  for (ni=0; ni<=ei; ni++)
  for (nj=0; nj<=ej; nj++)
  for (nk=0; nk<=ek; nk++)
  {  
    cell(l+1,2*i+ni,2*j+nj,2*k+nk)->setStatus(Leaf);
//    cell(l+1,2*i+ni,2*j+nj,2*k+nk)->setAverage(predictCell(l+1,2*i+ni,2*j+nj,2*k+nk));
  }

}  

/*
_______________________________________________________________________________

Combine cell
_______________________________________________________________________________

*/

void AdaptiveMesh::combine(const int l, const int i, const int j, const int k)
{
  int ni,nj,nk;
  
  int ei = 1;
  int ej = (Dimension > 1)? 1:0;
  int ek = (Dimension > 2)? 1:0;
  
  cell(l,i,j,k)->setStatus(Leaf);
  
  for (ni=0; ni<=ei; ni++)
  for (nj=0; nj<=ej; nj++)
  for (nk=0; nk<=ek; nk++)
    cell(l+1,2*i+ni,2*j+nj,2*k+nk)->setStatus(None);  

}  

/*
_______________________________________________________________________________

Predict value in cell
_______________________________________________________________________________

*/

Vector AdaptiveMesh::predictCell(const int l, const int i, const int j, const int k) const
{
	// --- Local variables ---

	int pi, pj=1, pk=1;	// Parity of Ni, Nj, Nk
	real a;

	Vector	result(QuantityNb);
	Vector V(QuantityNb);

	// --- Init result with the cell-average value of the father ---

	result = parent(l,i,j,k)->average();

	// --- 1D case ---

	pi = (i%2 == 0)?1:-1;
	a = pi*-.125;
	V += uncle(l,i,j,k,1,0,0)->average();
	V -= uncle(l,i,j,k,-1,0,0)->average();
	result += a*V;

	// --- 2D case ---

 	if (Dimension > 1)
	{
		pj = (j%2 == 0)?1:-1;
		a = pj*-.125;
		V.reset();
		V += uncle(l,i,j,k,0,1,0)->average();
		V -= uncle(l,i,j,k,0,-1,0)->average();
		result += a*V;
		
		a = pi*pj*.015625;
		V.reset();
		V += uncle(l,i,j,k,1,1,0)->average();
		V -= uncle(l,i,j,k,1,-1,0)->average();
		V -= uncle(l,i,j,k,-1,1,0)->average();
		V += uncle(l,i,j,k,-1,-1,0)->average();
		result += a*V;
	}

	// --- 3D case ---

 	if (Dimension > 2)
	{
		pk = (k%2 == 0)?1:-1;
		a = pk*-.125;
		V.reset();
		V += uncle(l,i,j,k,0,0,1)->average();
		V -= uncle(l,i,j,k,0,0,-1)->average();
		result += a*V;
		
		a = pi*pk*.015625;
		V.reset();
		V += uncle(l,i,j,k,1,0,1)->average();
		V -= uncle(l,i,j,k,1,0,-1)->average();
		V -= uncle(l,i,j,k,-1,0,1)->average();
		V += uncle(l,i,j,k,-1,0,-1)->average();
		result += a*V;
		
		a = pj*pk*.015625;
		V.reset();
		V += uncle(l,i,j,k,0,1,1)->average();
		V -= uncle(l,i,j,k,0,1,-1)->average();
		V -= uncle(l,i,j,k,0,-1,1)->average();
		V += uncle(l,i,j,k,0,-1,-1)->average();
		result += a*V;

		a = pi*pj*pk*-.001953125;
		V.reset();
		V += uncle(l,i,j,k,1,1,1)->average();
		V -= uncle(l,i,j,k,1,1,-1)->average();
		V -= uncle(l,i,j,k,1,-1,1)->average();
		V += uncle(l,i,j,k,1,-1,-1)->average();
		V -= uncle(l,i,j,k,-1,1,1)->average();
		V += uncle(l,i,j,k,-1,1,-1)->average();
		V += uncle(l,i,j,k,-1,-1,1)->average();
		V -= uncle(l,i,j,k,-1,-1,-1)->average();
		result += a*V;
	}
	return result;
}

/*
_______________________________________________________________________________

Project value in cell
_______________________________________________________________________________

*/

void AdaptiveMesh::projectCell(const int l, const int i, const int j, const int k)
{
  int ni, nj, nk;
  int ei = 1;
  int ej = (Dimension > 1)? 1:0;
  int ek = (Dimension > 2)? 1:0;
  real invN = 1./(1<<Dimension);
  
  cell(l,i,j,k)->resetAverage();
  
  for (ni=0; ni<=ei; ni++)
  for (nj=0; nj<=ej; nj++)
  for (nk=0; nk<=ek; nk++)
    cell(l,i,j,k)->setAverage(cell(l,i,j,k)->average()+cell(l+1,2*i+ni,2*j+nj,2*k+nk)->average());
  
  cell(l,i,j,k)->setAverage(invN*cell(l,i,j,k)->average());
}

/*
_______________________________________________________________________________

Update tree structure
_______________________________________________________________________________

*/

void AdaptiveMesh::update()
{
  int l,i,j,k;
  
  int ei=1;
  int ej=(Dimension > 1)? 1:0;
  int ek=(Dimension > 2)? 1:0;
  
  for (l=ScaleNb-1; l>=ScaleMin; l--)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=0; i<=ei*(N(l)-1); i++)
  for (j=0; j<=ej*(N(l)-1); j++)
  for (k=0; k<=ek*(N(l)-1); k++)
  {     
    if (cell(l,i,j,k)->status()==Node)
      if (isCousinOfLeaf(l,i,j,k) || isParentOfLeaf(l,i,j,k)) projectCell(l,i,j,k);
  }  
  
  for (l=ScaleMin+1; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=0; i<=ei*(N(l)-1); i++)
  for (j=0; j<=ej*(N(l)-1); j++)
  for (k=0; k<=ek*(N(l)-1); k++)
  { 
    if (cell(l,i,j,k)->status()==None)
      if (isCousinOfLeaf(l,i,j,k)) cell(l,i,j,k)->setAverage(predictCell(l,i,j,k));
  }  
} 

/*
_______________________________________________________________________________

Update complete tree structure 
_______________________________________________________________________________

*/

void AdaptiveMesh::updateAll()
{
  int l,i,j,k;
  
  int ei=1;
  int ej=(Dimension > 1)? 1:0;
  int ek=(Dimension > 2)? 1:0;
  
  for (l=ScaleNb-1; l>=ScaleMin; l--)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=0; i<=ei*(N(l)-1); i++)
  for (j=0; j<=ej*(N(l)-1); j++)
  for (k=0; k<=ek*(N(l)-1); k++)
  {     
    if (cell(l,i,j,k)->status()==Node)
      projectCell(l,i,j,k);
  }  
  
  for (l=ScaleMin+1; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=0; i<=ei*(N(l)-1); i++)
  for (j=0; j<=ej*(N(l)-1); j++)
  for (k=0; k<=ek*(N(l)-1); k++)
  { 
    if (cell(l,i,j,k)->status()==None)
      cell(l,i,j,k)->setAverage(predictCell(l,i,j,k));
  }  
} 

/*
_______________________________________________________________________________

Adapt tree structure
_______________________________________________________________________________

*/

void AdaptiveMesh::adapt()
{
  Vector D(QuantityNb);
  Vector d(QuantityNb);
  int smallChildNb;
  int totalChildNb = 1<<Dimension;
  int childLeafNb;
  
  int l,i,j,k;
  int ni, nj, nk;
  
  int ei=1;
  int ej=(Dimension > 1)? 1:0;
  int ek=(Dimension > 2)? 1:0;
  
  for (l=ScaleNb-1; l>=ScaleMin; l--)  
  for (i=0; i<=ei*(N(l)-1); i++)
  for (j=0; j<=ej*(N(l)-1); j++)
  for (k=0; k<=ek*(N(l)-1); k++)
  { 
    if (cell(l,i,j,k)->status()==Leaf)
      childLeafNb = totalChildNb;
    else
    {
      childLeafNb = 0;
      for (ni=0; ni<=ei; ni++)
      for (nj=0; nj<=ej; nj++)
      for (nk=0; nk<=ek; nk++)
	if (cell(l+1,2*i+ni,2*j+nj,2*k+nk)->status()==Leaf) childLeafNb++;
    }
    
    if (childLeafNb==totalChildNb)
    {
      // Compute details in the parent cell
    
      D = detail(l,i,j,k);
      smallChildNb = 0;
    
      // Compute details in the children cells
    
      for (ni=0; ni<=ei; ni++)
      for (nj=0; nj<=ej; nj++)
      for (nk=0; nk<=ek; nk++)
      {  
	d = detail(l+1,2*i+ni,2*j+nj,2*k+nk);
	if (isSmall(l+1,d)) smallChildNb++;
      }
    
      // If the detail in the local cell is big or if it is the uncle of a node, split
    
      if (cell(l,i,j,k)->status()==Leaf)
      {
	if (!isSmall(l,D)) split(l,i,j,k);
      }
      
      // If the detail is small in both parent and children, combine children
            
      if (isSmall(l,D) && smallChildNb==totalChildNb && !isUncleOfNode(l,i,j,k)) 
	combine(l,i,j,k);
    }
  }  

  for (l=ScaleNb-1; l>=ScaleMin; l--)  
  for (i=0; i<=ei*(N(l)-1); i++)
  for (j=0; j<=ej*(N(l)-1); j++)
  for (k=0; k<=ek*(N(l)-1); k++)
  { 
    if (cell(l,i,j,k)->status()==Leaf)
      if (isUncleOfNode(l,i,j,k)) split(l,i,j,k);
  }
}

/*
_______________________________________________________________________________

Compute leaves and cells
_______________________________________________________________________________

*/

void AdaptiveMesh::countLeavesAndCells()
{

  int l, i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  NodeNb=0;
  LeafNb=0;

  // Loop on real cells
  
  for (l=0; l<=ScaleNb; l++)
  for (i=0; i <= ei*(N(l) - 1); i++)
  for (j=0; j <= ej*(N(l) - 1); j++)
  for (k=0; k <= ek*(N(l) - 1); k++)
  {
    if (cell(l,i,j,k)->status()==Leaf) LeafNb++;
    if (cell(l,i,j,k)->status()==Leaf || cell(l,i,j,k)->status()==Node) NodeNb++;
  }
}

/*
_______________________________________________________________________________

Compute detail
_______________________________________________________________________________

*/

Vector AdaptiveMesh::detail(const int l, const int i, const int j, const int k) const
{
    Vector result(QuantityNb);
    Vector Q(QuantityNb);
    Vector Qp(QuantityNb);
    Vector M(Dimension);
    Vector Mp(Dimension);
    real NormM, NormMp;
    
    int axisNo;
    int Nl = 1<<ScaleNb;

    switch(Equation)
    { 
      case ConvectionDiffusion:
      default:
	return (cell(l,i,j,k)->average()-predictCell(l,i,j,k));
	
      case Burgers:
	
	// For the detail, take the norm of the velocity
	
	M = cell(l,i,j,k)->average();
	Mp = predictCell(l,i,j,k);
	NormM = sqrt(M*M);
	NormMp = sqrt(Mp*Mp);
	
	for (axisNo=1; axisNo<=Dimension; axisNo++)
	{  
	  M.set(axisNo,NormM);
	  Mp.set(axisNo,NormMp);
	}
	
	return M-Mp;

      case NavierStokes:
	Q = cell(l,i,j,k)->average();
	Qp = predictCell(l,i,j,k);
	
	// For the detail, take the norm of the momentum 
	
	
	for (axisNo=1; axisNo<=Dimension; axisNo++)
	{  
	  M.set(axisNo,Q.at(axisNo+1));
	  Mp.set(axisNo,Qp.at(axisNo+1));
	}
	
	NormM = sqrt(M*M);
	NormMp = sqrt(Mp*Mp);
	
	for (axisNo=1; axisNo<=Dimension; axisNo++)
	{  
	  Q.set(axisNo+1,NormM);
	  Qp.set(axisNo+1,NormMp);
	}
	
	// For the Sedov problem, refine in the center to avoid numerical instability
	
	if (RefineInCenter)
	{
	  switch (Dimension)
	  {
	    case 1:
	      if (l == ScaleNb && SQR(i-Nl/2) <= SQR(RefineWidth/2))
		  Qp.set(1,-1.E+99);
	      break;
	   
	    case 2:
	      if (l == ScaleNb && SQR(i-Nl/2)+SQR(j-Nl/2) <= SQR(RefineWidth/2))
		Qp.set(1,-1.E+99);
	      break;
	    
	    case 3:
	      if (l == ScaleNb && SQR(i-Nl/2)+SQR(j-Nl/2)+ SQR(k-Nl/2) <= SQR(RefineWidth/2))
		Qp.set(1,-1.E+99);
	      break;
	  };
	}
	
	
	return Q-Qp;
    };
}

/*
_______________________________________________________________________________

returns true if (l,i,j,k) is the uncle of a node (and needs to be split)
_______________________________________________________________________________

*/

bool AdaptiveMesh::isUncleOfNode(const int l, const int i, const int j, const int k) const
{
  int ni, nj, nk;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  if (l >= ScaleNb) return false;
  if (cell(l,i,j,k)->status() !=Leaf) return false;
  
  for (ni=-ei; ni<=ei; ni++)
  for (nj=-ej; nj<=ej; nj++)
  for (nk=-ek; nk<=ek; nk++)
  {  
    if (!(ni==0 && nj==0 && nk==0) && cell(l+1,2*(i+ni),2*(j+nj),2*(k+nk))->status()==Node)
      return true;
  }
  return false;
}

/*
_______________________________________________________________________________

returns true if (l,i,j,k) is the cousin of a leaf
_______________________________________________________________________________

*/

bool AdaptiveMesh::isCousinOfLeaf(const int l, const int i, const int j, const int k) const
{
  int ni, nj, nk;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  if (l >= ScaleNb) return false;
  
  for (ni=-3*ei; ni<=3*ei; ni++)
  for (nj=-3*ej; nj<=3*ej; nj++)
  for (nk=-3*ek; nk<=3*ek; nk++)
  {  
    if (cell(l,i+ni,j+nj,k+nk)->status()==Leaf)
      return true;
  }
  return false;
}

/*
_______________________________________________________________________________

returns true if (l,i,j,k) is the parent of a leaf
_______________________________________________________________________________

*/

bool AdaptiveMesh::isParentOfLeaf(const int l, const int i, const int j, const int k) const
{
  int ni, nj, nk;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  if (l >= ScaleNb) return false;
  
  for (ni=0; ni<=ei; ni++)
  for (nj=0; nj<=ej; nj++)
  for (nk=0; nk<=ek; nk++)
  {  
    if (cell(l+1,i+ni,j+nj,k+nk)->status()==Leaf)
      return true;
  }
  return false;
}
/*
_______________________________________________________________________________

returns true if the vector average is smaller than the tolerance
_______________________________________________________________________________

*/

bool AdaptiveMesh::isSmall(const int l, const Vector& average) const
{
    int i;
    real eps = exp((Dimension*(l-ScaleNb))*log(2.))*Tolerance;

//    eps *= (IterationNo < 1)? 1.: (1.-IterationNo*exp(-1.))/IterationNo;
    
    for (i=1; i<=QuantityNb; i++)
      if (fabs(average.at(i)/MaxAverage.at(i)) >= eps) return false;
      
    return true;
} 

/*
_______________________________________________________________________________

Copy cell-average values from sourceStageNo to targetStageNo
_______________________________________________________________________________

*/

void AdaptiveMesh::copyAverage(const int sourceStorageNo, const int targetStorageNo)
{
  int l, i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  // Loop on real cells
  
  for (l=ScaleMin; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=-3*ei; i <= ei*(N(l) + 2); i++)
  for (j=-3*ej; j <= ej*(N(l) + 2); j++)
  for (k=-3*ek; k <= ek*(N(l) + 2); k++)
  {
    if (cell(l,i,j,k)->status()==Leaf)
      cell(l,i,j,k)->copyStorage(sourceStorageNo, targetStorageNo);
  }
} 


/*
_______________________________________________________________________________

Set boundary conditions
_______________________________________________________________________________

*/

void AdaptiveMesh::setBoundaryConditions()
{
  int l, i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  int Nl;
  int imin,jmin,kmin;
  int imax,jmax,kmax;
  
  i=j=k=0;
  imin=jmin=kmin=0;
  
  // Boundary conditions in XMin
  
  for (l=ScaleMin; l<=ScaleNb; l++)
  {
    Nl = N(l);
    imax = ei*(Nl-1);
    jmax = ej*(Nl-1);
    kmax = ek*(Nl-1);
    
    // --- XMin ---------------------------------------------------------------

#pragma omp parallel for private(j,k) shared(l,Nl,jmin,jmax,kmin,kmax)
    for (j=jmin; j <= jmax; j++)
    for (k=kmin; k <= kmax; k++)
    {
	cell(l,-3,j,k)->linkBoundaryCell(BCXMin,cell(l,2,j,k),cell(l,-2,j,k),cell(l,Nl-3,j,k),1);
	cell(l,-2,j,k)->linkBoundaryCell(BCXMin,cell(l,1,j,k),cell(l,-1,j,k),cell(l,Nl-2,j,k),1);
	cell(l,-1,j,k)->linkBoundaryCell(BCXMin,cell(l,0,j,k),cell(l, 0,j,k),cell(l,Nl-1,j,k),1);
    }
    
    // --- XMax ---------------------------------------------------------------
        
#pragma omp parallel for private(j,k) shared(l,Nl,jmin,jmax,kmin,kmax)
    for (j=jmin; j <= jmax; j++)
    for (k=kmin; k <= kmax; k++)
    {
	  cell(l,Nl  ,j,k)->linkBoundaryCell(BCXMax,cell(l,Nl-1,j,k),cell(l,Nl-1,j,k),cell(l,0,j,k),1);
	  cell(l,Nl+1,j,k)->linkBoundaryCell(BCXMax,cell(l,Nl-2,j,k),cell(l,Nl  ,j,k),cell(l,1,j,k),1);
	  cell(l,Nl+2,j,k)->linkBoundaryCell(BCXMax,cell(l,Nl-3,j,k),cell(l,Nl+1,j,k),cell(l,2,j,k),1);
    }
 
  
    if (Dimension > 1)
    {
      // --- YMin -------------------------------------------------------------
      
#pragma omp parallel for private(i,k) shared(l,Nl,imin,imax,kmin,kmax)
	for (i=imin; i <= imax; i++)
	for (k=kmin; k <= kmax; k++)
	{
	    cell(l,i,-3,k)->linkBoundaryCell(BCYMin,cell(l,i,2,k),cell(l,i,-2,k),cell(l,i,Nl-3,k),2);
	    cell(l,i,-2,k)->linkBoundaryCell(BCYMin,cell(l,i,1,k),cell(l,i,-1,k),cell(l,i,Nl-2,k),2);
	    cell(l,i,-1,k)->linkBoundaryCell(BCYMin,cell(l,i,0,k),cell(l,i, 0,k),cell(l,i,Nl-1,k),2);
	}
    
      
      // --- YMax -------------------------------------------------------------
      
#pragma omp parallel for private(i,k) shared(l,Nl,imin,imax,kmin,kmax)
	for (i=imin; i <= imax; i++)
	for (k=kmin; k <= kmax; k++)
	{
	    cell(l,i,Nl  ,k)->linkBoundaryCell(BCYMax,cell(l,i,Nl-1,k),cell(l,i,Nl-1,k),cell(l,i,0,k),2);
	    cell(l,i,Nl+1,k)->linkBoundaryCell(BCYMax,cell(l,i,Nl-2,k),cell(l,i,Nl  ,k),cell(l,i,1,k),2);
	    cell(l,i,Nl+2,k)->linkBoundaryCell(BCYMax,cell(l,i,Nl-3,k),cell(l,i,Nl+1,k),cell(l,i,2,k),2);
	}
      
    }
   
    if (Dimension > 2)
    {
      
      // --- ZMin -------------------------------------------------------------
      
#pragma omp parallel for private(i,j) shared(l,Nl,imin,imax,jmin,jmax)
	for (i=imin; i <= imax; i++)
	for (j=jmin; j <= jmax; j++)
	{
	    cell(l,i,j,-3)->linkBoundaryCell(BCZMin,cell(l,i,j,2),cell(l,i,j,-2),cell(l,i,j,Nl-3),3);
	    cell(l,i,j,-2)->linkBoundaryCell(BCZMin,cell(l,i,j,1),cell(l,i,j,-1),cell(l,i,j,Nl-2),3);
	    cell(l,i,j,-1)->linkBoundaryCell(BCZMin,cell(l,i,j,0),cell(l,i,j, 0),cell(l,i,j,Nl-1),3);
	}
      
      // --- ZMax -------------------------------------------------------------
      
#pragma omp parallel for private(i,j) shared(l,Nl,imin,imax,jmin,jmax)
	for (i=imin; i <= imax; i++)
	for (j=jmin; j <= jmax; j++)
	{
	    cell(l,i,j,Nl  )->linkBoundaryCell(BCZMax,cell(l,i,j,Nl-1),cell(l,i,j,Nl-1),cell(l,i,j,0),3);
	    cell(l,i,j,Nl+1)->linkBoundaryCell(BCZMax,cell(l,i,j,Nl-2),cell(l,i,j,Nl  ),cell(l,i,j,1),3);
	    cell(l,i,j,Nl+2)->linkBoundaryCell(BCZMax,cell(l,i,j,Nl-3),cell(l,i,j,Nl+1),cell(l,i,j,2),3);
	}
      
    }    
  }
}

/*
_______________________________________________________________________________

compute divergence
_______________________________________________________________________________

*/

void AdaptiveMesh::computeDivergence()
{
  int l, i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  Vector Fx1,Fx2,Fy1,Fy2,Fz1,Fz2;
  
  // Loop on cells

  for (l=ScaleMin; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k,Fx1,Fx2,Fy1,Fy2,Fz1,Fz2) shared(l,ei,ej,ek)
  for (i=0; i <= ei*(N(l)-1); i++)
  for (j=0; j <= ej*(N(l)-1); j++)
  for (k=0; k <= ek*(N(l)-1); k++)
  {    
    if (cell(l,i,j,k)->status()==Leaf)
    {            
      // Compute fluxes in the x-direction
    
      Fx2 = flux(cell(l,i-2,j,k),cell(l,i-1,j,k),cell(l,i,j,k),cell(l,i+1,j,k),cell(l,i+2,j,k),cell(l,i+3,j,k),1);
      Fx1 = flux(cell(l,i-3,j,k),cell(l,i-2,j,k),cell(l,i-1,j,k),cell(l,i,j,k),cell(l,i+1,j,k),cell(l,i+2,j,k),1);
    
      cell(l,i,j,k)->setDivergence((Fx1-Fx2)/cell(l,i,j,k)->dx());
    
      if (Dimension > 1)
      {
	// Compute fluxes in the y-direction
      
	Fy2 = flux(cell(l,i,j-2,k),cell(l,i,j-1,k),cell(l,i,j,k),cell(l,i,j+1,k),cell(l,i,j+2,k),cell(l,i,j+3,k),2);
	Fy1 = flux(cell(l,i,j-3,k),cell(l,i,j-2,k),cell(l,i,j-1,k),cell(l,i,j,k),cell(l,i,j+1,k),cell(l,i,j+2,k),2);
    
      cell(l,i,j,k)->addDivergence((Fy1-Fy2)/cell(l,i,j,k)->dy());   
      }

      if (Dimension > 2)
      {
	// Compute fluxes in the z-direction
      
	Fz2 = flux(cell(l,i,j,k-2),cell(l,i,j,k-1),cell(l,i,j,k),cell(l,i,j,k+1),cell(l,i,j,k+2),cell(l,i,j,k+3),3);
	Fz1 = flux(cell(l,i,j,k-3),cell(l,i,j,k-2),cell(l,i,j,k-1),cell(l,i,j,k),cell(l,i,j,k+1),cell(l,i,j,k+2),3);
    
	cell(l,i,j,k)->addDivergence((Fz1-Fz2)/cell(l,i,j,k)->dz());   
      }
    }
  }
}

/*
_______________________________________________________________________________

compute Runge-Kutta step
_______________________________________________________________________________

*/

void AdaptiveMesh::computeStep()
{
  int l, i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  real X,Y,Z,t;  
  
  const real a = 1./3.;

  Vector U0, U, D;
  real dt = TimeStep;

  X=Y=Z=t=0.;
  
  // Loop on cells
  
  for (l=ScaleMin; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k,U0,U,D,X,Y,Z,t) shared(l,ei,ej,ek,dt)
  for (i=0; i <= ei*(N(l)-1); i++)
  for (j=0; j <= ej*(N(l)-1); j++)
  for (k=0; k <= ek*(N(l)-1); k++)
  {
    if (cell(l,i,j,k)->status()==Leaf)
    {
      
      switch (StepNo)
      {
	case 1:
	default:	
	  U = cell(l,i,j,k)->average();
	  D = cell(l,i,j,k)->divergence();
	  cell(l,i,j,k)->setAverage(U + dt*D);
	  break;

	case 2:
	  switch(StepNb)
	  {
	    case 1:
	    default:
	      break;
      
	    case 2:
	      U0 = cell(l,i,j,k)->average(2);
	      U  = cell(l,i,j,k)->average();
	      D  = cell(l,i,j,k)->divergence();
	      cell(l,i,j,k)->setAverage(0.5*(U0 + U + dt*D) );
	      break;
	    
	    case 3:
	      U0 = cell(l,i,j,k)->average(2);
	      U  = cell(l,i,j,k)->average();
	      D  = cell(l,i,j,k)->divergence();
	      cell(l,i,j,k)->setAverage(0.25*(3.*U0 + U + dt*D) );
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
	      U0 = cell(l,i,j,k)->average(2);
	      U  = cell(l,i,j,k)->average();
	      D  = cell(l,i,j,k)->divergence();
	      cell(l,i,j,k)->setAverage(a*(U0 + 2.*(U + dt*D)) );
	      break;    
	  };
	  break;
      }; 
      
      if (UsePenalization)
      {
	X = cell(l,i,j,k)->x();
	Y = cell(l,i,j,k)->y();
	Z = cell(l,i,j,k)->z();
	t = TotalTime-RemainingTime;
	cell(l,i,j,k)->setAverage(QuantityNb,1,(isImmersed(X,Y,Z,t)? 0.:1.));
      }
      
    }
  }
}
/*
_______________________________________________________________________________

Compute integral values
_______________________________________________________________________________

*/

void AdaptiveMesh::computeIntegralValues()
{
  int l,i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  int axisNo,speciesNo;
  real dx, dy, dz;
  real x,y,z;
  int ll,li,lj,lk;
  real pred;	// prediction error
  real mu;	// dynamic viscosity
    
  MaxEigenvalue = 0.;
  KineticEnergy = 0.;
  MaxDiffusivity = (Equation==NavierStokes)? 0.: Diffusivity;
  MaxAverage.reset();
  PredictionError.reset();

  if (IterationNo > 0)
  {
    AverageTimeStep = ((IterationNo-1)*AverageTimeStep + TimeStep)/IterationNo;
    AverageNodeNb = ((IterationNo-1)*AverageNodeNb + NodeNb)/IterationNo;
    AverageLeafNb = ((IterationNo-1)*AverageLeafNb + LeafNb)/IterationNo;
  }
  else
  {
    AverageTimeStep = TimeStep;
    AverageNodeNb = NodeNb;
    AverageLeafNb = LeafNb;
  }
  // Loop on real cells
  
  for (l=ScaleMin; l<=ScaleNb; l++)
  for (i=0; i <= ei*(N(l)-1); i++)
  for (j=0; j <= ej*(N(l)-1); j++)
  for (k=0; k <= ek*(N(l)-1); k++)
  {  
    if (l==ScaleNb)
    {
 
      ll=ScaleNb-1;
      li=i/2;
      lj=j/2;
      lk=k/2;
       
      // Compute difference between finest grid solution and next coarser grid solution
      
      for (axisNo=1;axisNo<=QuantityNb;axisNo++)
	PredictionError.set(axisNo,MAX(PredictionError.at(axisNo),ABS(cell(l,i,j,k)->average(axisNo,1)-cell(ll,li,lj,lk)->average(axisNo,1))));
 //	PredictionError.set(axisNo,PredictionError.at(axisNo)+ABS(cell(l,i,j,k)->average(axisNo,1)-cell(ll,li,lj,lk)->average(axisNo,1)));
   }
    
    if (cell(l,i,j,k)->status()==Leaf)
    {  
      MaxEigenvalue = MAX(MaxEigenvalue,cell(l,i,j,k)->maxEigenvalue());
    
      if (Equation == NavierStokes)
      {
	mu = cell(l,i,j,k)->sutherland(cell(l,i,j,k)->T());
	MaxDiffusivity = MAX(MaxDiffusivity, MAX(Pr,1.)*mu/cell(l,i,j,k)->rho());    
      }  
    
      dx = cell(l,i,j,k)->dx();
      dy = (Dimension > 1)? cell(l,i,j,k)->dy():1.;
      dz = (Dimension > 2)? cell(l,i,j,k)->dz():1.;
      x = cell(l,i,j,k)->x();
      y = (Dimension > 1) ? cell(l,i,j,k)->y():0.;
      z = (Dimension > 2) ? cell(l,i,j,k)->z():0.;

      // Compute depending on the equation
      
      switch(Equation)
      { 
	case ConvectionDiffusion:
	  KineticEnergy += ABS(cell(l,i,j,k)->density()-initCondition(x,y,z).at(1))*dx*dy*dz;
	  MaxAverage.set(1,MAX(MaxAverage.at(1),ABS(cell(l,i,j,k)->density())));
	  break;
	
	case Burgers:
	  KineticEnergy += 0.5*SQR(cell(l,i,j,k)->velocityNorm())*dx*dy*dz;
	  for (axisNo=1; axisNo<=Dimension; axisNo++)
	    MaxAverage.set(axisNo,MAX(MaxAverage.at(axisNo),cell(l,i,j,k)->velocityNorm()));
	  break;
	
	case NavierStokes:
	  KineticEnergy += 0.5*cell(l,i,j,k)->density()*SQR(cell(l,i,j,k)->velocityNorm())*dx*dy*dz;
	  MaxAverage.set(1,MAX(MaxAverage.at(1),cell(l,i,j,k)->density()));
	  for (axisNo=1; axisNo<=Dimension; axisNo++)
	    MaxAverage.set(axisNo+1,MAX(MaxAverage.at(axisNo+1),cell(l,i,j,k)->momentumNorm()));
	  MaxAverage.set(Dimension+2,MAX(MaxAverage.at(1),ABS(cell(l,i,j,k)->energy())));
	  if (SpeciesNb > 1)
	  {
	    for (speciesNo = 1; speciesNo <= SpeciesNb-1; speciesNo++)
	      MaxAverage.set(Dimension+2+speciesNo,MAX(MaxAverage.at(Dimension+2+speciesNo),cell(l,i,j,k)->average(Dimension+2+speciesNo,1)));
	  }
	  
	  break;
      };
    }
  }
  
  // Compute prediction error
  
  pred = 0.;
  
  for (axisNo=1;axisNo<=QuantityNb;axisNo++)
    pred = MAX(pred, PredictionError.at(axisNo)/(Tolerance*MaxAverage.at(axisNo)));
//    pred = MAX(pred, PredictionError.at(axisNo)/((1<<(ScaleNb*Dimension))*Tolerance*MaxAverage.at(axisNo)));
  
  PredictionError.set(1,pred);
  
  // Compute the maximum of the density averaged in time (for detonation only)
  
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

void AdaptiveMesh::computeTimeStep()
{
  real K = MultiresolutionError;
    
  if (CFL != 0.)
  {

    if (IsInviscid)
      TimeStep = CFL*SpaceStep/MaxEigenvalue;
    else
      TimeStep = CFL*MIN(SpaceStep/MaxEigenvalue, SQR(SpaceStep)/(2.*MaxDiffusivity));
  
    if (TimeStep > RemainingTime) TimeStep = RemainingTime;
  }

  if (UseVariableTolerance)
  {
    if (IterationNo > 1)
      Tolerance *= (1.-InitialTolerance/K);
  
  }
  
}

/*
_______________________________________________________________________________

Compute gradient
_______________________________________________________________________________

*/

void AdaptiveMesh::computeGradient()
{
  int l, i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  // Loop on real cells
  
  for (l=ScaleMin; l<=ScaleNb;l++)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=0; i <= ei*(N(l)-1); i++)
  for (j=0; j <= ej*(N(l)-1); j++)
  for (k=0; k <= ek*(N(l)-1); k++)
  {
    if (cell(l,i,j,k)->status()==Leaf)
    {  
      cell(l,i,j,k)->computeGradient(cell(l,i-2,j,k),cell(l,i-1,j,k),cell(l,i+1,j,k),cell(l,i+2,j,k),1);
      if (Dimension > 1) cell(l,i,j,k)->computeGradient(cell(l,i,j-2,k),cell(l,i,j-1,k),cell(l,i,j+1,k),cell(l,i,j+2,k),2);
      if (Dimension > 2) cell(l,i,j,k)->computeGradient(cell(l,i,j,k-2),cell(l,i,j,k-1),cell(l,i,j,k+1),cell(l,i,j,k+2),3);
    }
  }
}

/*
_______________________________________________________________________________

Compute source
_______________________________________________________________________________

*/

void AdaptiveMesh::computeSource()
{
  int l, i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  // Loop on real cells
  
  for (l=ScaleMin; l<=ScaleNb;l++)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=0; i <= ei*(N(l)-1); i++)
  for (j=0; j <= ej*(N(l)-1); j++)
  for (k=0; k <= ek*(N(l)-1); k++)
  {
    if (cell(l,i,j,k)->status()==Leaf)
      cell(l,i,j,k)->computeSource();
  }
}

/*
_______________________________________________________________________________

Compute initial condition
_______________________________________________________________________________

*/

void AdaptiveMesh::initAverage()
{
  int l,i,j,k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  real x,y,z;
  
  // Fill tree
  
  for (l=0; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k,x,y,z) shared(l,ei,ej,ek)
  for (i=-3*ei; i<=ei*(N(l)+2); i++)
  for (j=-3*ej; j<=ej*(N(l)+2); j++)
  for (k=-3*ek; k<=ek*(N(l)+2); k++)
  {
    x = cell(l,i,j,k)->x();
    y = (Dimension > 1) ? cell(l,i,j,k)->y() : 0.;
    z = (Dimension > 2) ? cell(l,i,j,k)->z() : 0.;
    cell(l,i,j,k)->setAverage(1,initCondition(x,y,z));

    if (i<0 || i>=N(l) || j<0 || j>=N(l) || k<0 || k>=N(l))
      cell(l,i,j,k)->setStatus(Boundary);
    else if (l==ScaleNb)
      cell(l,i,j,k)->setStatus(Leaf);
    else
      cell(l,i,j,k)->setStatus(Node);
  }  
}

/*
_______________________________________________________________________________

Check the stability of the computation
_______________________________________________________________________________

*/

void AdaptiveMesh::checkStability()
{
  int l, i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  FILE* f;
  
  // Loop on real cells
  
  for (l = ScaleMin; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k,f) shared(l,ei,ej,ek)
  for (i=0; i <= ei*(N(l)-1); i++)
  for (j=0; j <= ej*(N(l)-1); j++)
  for (k=0; k <= ek*(N(l)-1); k++)
  {
    if (cell(l,i,j,k)->status()==Leaf)
    {  
      if (cell(l,i,j,k)->isNaN())
      {
	printf("lara: Numerical instability detected at iteration %i, position (%15.8f,%15.8f,%15.8f)\n", IterationNo, cell(l,i,j,k)->x(), (Dimension > 1) ? cell(l,i,j,k)->y():0.,
	     (Dimension > 2) ? cell(l,i,j,k)->z():0.);
	f=fopen("performance.out","a");
	fprintf(f,"Numerical instability detected at iteration %i, position (%15.8f,%15.8f,%15.8f)\n\n", IterationNo, cell(l,i,j,k)->x(), (Dimension > 1) ? cell(l,i,j,k)->y():0.,
	     (Dimension > 2) ? cell(l,i,j,k)->z():0.);
	fclose(f);
	exit(-1);
      }
    }  
  }
}

/*
_______________________________________________________________________________

Write cell-average values
_______________________________________________________________________________

*/

void AdaptiveMesh::writeAverage()
{
  int Nl = N(ScaleNb);
  int l = ScaleNb;

  int i,j,k,n;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  int imin = 0;
  int imax = ei*(Nl-1);
  int jmin = 0;
  int jmax = ej*(Nl-1);
  int kmin = 0;
  int kmax = ek*(Nl-1);
  
  real dx = (XMax-XMin)/Nl;
  real dy = (YMax-YMin)/Nl;
  real dz = (ZMax-ZMin)/Nl; 
  
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
      fprintf(f,"DIMENSIONS %i %i %i\n", Nl+1, (Dimension < 2)? 2:Nl+1, (Dimension < 3)? 2:Nl+1);
      fprintf(f,"ORIGIN %f %f %f\n", XMin*ei, YMin*ej, ZMin*ek);
      fprintf(f,"SPACING %f %f %f\n", dx, (Dimension < 2)? 1.:dy, (Dimension < 3) ? 1.:dz);
      fprintf(f,"CELL_DATA %i\n", (int)(pow(Nl,Dimension)));

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
	  writeLn(f,cell(l,i,j,k)->rho());
      }    
      
      if (Equation == NavierStokes)
      {
	if (WriteType==Binary) fprintf(f,"\n");
	fprintf(f,"SCALARS p float\n");
	fprintf(f,"LOOKUP_TABLE default\n");
      
	for (k=kmin; k<=kmax; k++)
	for (j=jmin; j<=jmax; j++)
	for (i=imin; i<=imax; i++)
	  writeLn(f,cell(l,i,j,k)->p());
      
	if (WriteType==Binary) fprintf(f,"\n");
	fprintf(f,"SCALARS T float\n");
	fprintf(f,"LOOKUP_TABLE default\n");
      
	for (k=kmin; k<=kmax; k++)
	for (j=jmin; j<=jmax; j++)
	for (i=imin; i<=imax; i++)
	  writeLn(f,cell(l,i,j,k)->T());
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
	    writeLn(f,cell(l,i,j,k)->v(1));
	}
	else
	{
	  if (WriteType==Binary) fprintf(f,"\n");
	  fprintf(f,"VECTORS v float\n");
	
	  for (k=kmin; k<=kmax; k++)
	  for (j=jmin; j<=jmax; j++)
	  for (i=imin; i<=imax; i++)
	  {
	    write(f,cell(l,i,j,k)->v(1));
	    write(f,cell(l,i,j,k)->v(2));
	    write(f,(Dimension > 2) ? cell(l,i,j,k)->v(3):0.);
	    if (WriteType==ASCII) fprintf(f,"\n"); 
	  }
	}
      }
      
      if (Equation == NavierStokes && SpeciesNb > 1)
      {
	for (n=1; n<=SpeciesNb-1; n++)
	{
	  if (WriteType==Binary) fprintf(f,"\n");
	  
	  if (UsePenalization && n==SpeciesNb-1)
	    fprintf(f,"SCALARS Mask float\n");
	  else
	    fprintf(f,"SCALARS Y%i float\n",n);
	  
	  fprintf(f,"LOOKUP_TABLE default\n");
	  
	  for (k=kmin; k<=kmax; k++)
	  for (j=jmin; j<=jmax; j++)
	  for (i=imin; i<=imax; i++)
	  {
	    if (UsePenalization && n==SpeciesNb-1)
	      writeLn(f,cell(l,i,j,k)->average(QuantityNb,1));
	    else
	      writeLn(f,cell(l,i,j,k)->Y(n)); 
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
      for (n=1; n<=SpeciesNb-1; n++)
	fprintf(f,"             Y%i ",n); 
    }
    fprintf(f,"\n");

    // --- DATA ---
    
    for (k=kmin; k<=kmax; k++)
    {
      for (j=jmin; j<=jmax; j++)
      {
	for (i=imin; i<=imax; i++)
	{
	  for (n=1; n<=Dimension;n++)
	    write(f,cell(l,i,j,k)->center(n));
	      
	  if (Equation != Burgers)    
	    write(f,cell(l,i,j,k)->rho());
	  
	  if (Equation == NavierStokes)
	  {
	    write(f,cell(l,i,j,k)->p());
	    write(f,cell(l,i,j,k)->T());
	  }
	  
	  if (Equation != ConvectionDiffusion)
	  {
	    for (n=1; n<=Dimension;n++)
	      write(f,cell(l,i,j,k)->v(n));
	  }

	  if (Equation == NavierStokes && SpeciesNb > 1)
	  {
	    for (n=1; n<=SpeciesNb-1;n++)
	      write(f,cell(l,i,j,k)->Y(n));
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

Write cell-average values into columns (for error analysis)
_______________________________________________________________________________

*/

void AdaptiveMesh::writeStrip()
{
  int Nl = N(ScaleNb);
  int i,j,k,l;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  int imin = 0;
  int imax = ei*(Nl-1);
  int jmin = 0;
  int jmax = ej*(Nl-1);
  int kmin = 0;
  int kmax = ek*(Nl-1);
  
  char fileName[256];
  FILE *f;
  
  // Write data
  
  sprintf(fileName,"Strip%05i.gnu",WriteNo-1);
  
  // Open file
  
  f = fopen(fileName,"w");
    
  for (k=kmin; k<=kmax; k++)
  for (j=jmin; j<=jmax; j++) 
  for (i=imin; i<=imax; i++)
  {
      for (l=1; l<= QuantityNb; l++)
	write(f,cell(ScaleNb,i,j,k)->average(l,1));
      fprintf(f,"\n");
  }

  
  fclose(f);
  
}
/*
_______________________________________________________________________________

Write tree structure
_______________________________________________________________________________

*/

void AdaptiveMesh::writeTree()
{
  int l,i,j,k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  int ni, nj, nk;
  int index;
    
  char fileName[256];
  FILE *f;
  
  // Write data
  
  if (WriteFormat == VTK)
    sprintf(fileName,"Tree%05i.vtk",WriteNo);
  else
    sprintf(fileName,"Tree%05i.gnu",WriteNo);
  
  // Open file
  
  f = fopen(fileName,"w");
    
  switch (WriteFormat)
  {  
    case VTK:
      
      // --- HEADER ---
      
      fprintf(f,"# vtk DataFile Version 3.0\n");
      fprintf(f,"Generated by Lara\n");
      fprintf(f,"ASCII\n");
      fprintf(f,"DATASET POLYDATA\n");
      
      // --- Write points ---
      
      fprintf(f,"POINTS %i float\n", (1<<Dimension)*LeafNb);
      for (l=ScaleMin; l<=ScaleNb;l++)
      for (k=0; k<=ek*(N(l)-1); k++)
      for (j=0; j<=ej*(N(l)-1); j++)
      for (i=0; i<=ei*(N(l)-1); i++)
      {
	if (cell(l,i,j,k)->status()==Leaf)
	{
	  for (nk=0; nk<=ek; nk++)
	  for (nj=0; nj<=ej; nj++)
	  for (ni=0; ni<=ei; ni++)
	  {
	    fprintf(f,"%15.8e ",cell(l,i,j,k)->x()+(ni-0.5)*cell(l,i,j,k)->dx());
	    fprintf(f,"%15.8e ",(Dimension > 1) ? cell(l,i,j,k)->y()+(nj-0.5)*cell(l,i,j,k)->dy():0.);
	    fprintf(f,"%15.8e ",(Dimension > 2) ? cell(l,i,j,k)->z()+(nj-0.5)*cell(l,i,j,k)->dz():0.);
	    fprintf(f,"\n");
	  }
	}
      }

      // --- Write cells ---
      
      if (Dimension == 1)
	fprintf(f,"LINES %i %i\n",LeafNb,3*LeafNb);
      else if (Dimension == 2)
	fprintf(f,"POLYGONS %i %i\n",LeafNb,5*LeafNb);
      else
	fprintf(f,"POLYGONS %i %i\n",4*LeafNb,20*LeafNb);
      
      index = 0;
      for (i=0;i<=LeafNb-1; i++)
      {
	switch(Dimension)
	{
	  case 1:
	  default:
	    fprintf(f,"2 %15i %15i\n", index, index+1);
	    break;
	    
	  case 2:
	    fprintf(f,"4 %15i %15i %15i %15i\n", index, index+1, index+3, index+2);
	    break;
	    
	  case 3:
	    fprintf(f,"4 %15i %15i %15i %15i\n", index  , index+1, index+3, index+2);
	    fprintf(f,"4 %15i %15i %15i %15i\n", index+4, index+5, index+7, index+6);
	    fprintf(f,"4 %15i %15i %15i %15i\n", index  , index+1, index+5, index+4);
	    fprintf(f,"4 %15i %15i %15i %15i\n", index+2, index+3, index+7, index+6);
	    break;
	}

	index += (1<<Dimension);
	fprintf(f,"\n");
      }
      
      break;
    
    case Gnuplot:
   
      // --- HEADER ---
    
      fprintf(f, "# %13s ", "x");
      if (Dimension > 1) fprintf(f, "%15s ", "y");
      if (Dimension > 2) fprintf(f, "%15s ", "z");
      fprintf(f, "%15s\n", "l");
    
      // --- DATA ---
    
      for (l=ScaleMin; l<=ScaleNb;l++)
      for (k=0; k<=ek*(N(l)-1); k++)
      for (j=0; j<=ej*(N(l)-1); j++)
      for (i=0; i<=ei*(N(l)-1); i++)
      {
	if (cell(l,i,j,k)->status()==Leaf)
	{
	  fprintf(f,"%15.8e ",cell(l,i,j,k)->x());
	  if (Dimension > 1) fprintf(f,"%15.8e ",cell(l,i,j,k)->y());
	  if (Dimension > 2) fprintf(f,"%15.8e ",cell(l,i,j,k)->z());
	  fprintf(f,"%15i\n",l);
	}
      }
      break;  
  };
	  
  fclose(f);
  
}

/*
_______________________________________________________________________________

Write integral values
_______________________________________________________________________________

*/

void AdaptiveMesh::writeIntegral()
{
  FILE *f;
  
  real delta = sqrt(2*(Gamma-1.)/(Gamma+1.));
  real rhoVN = Gamma/(1-delta);
  real err = ABS(rhoVN-TimeAverageMaxDensity)/ABS(rhoVN);
  
  if (IterationNo == 0 && !Restart)
  {  
    f = fopen("Integral.gnu","w");
    fprintf(f,"# %13s %15s %15s %15s %15s %15s %15s %15s\n", "Time", "Time step", "Kin. energy", "Max density", "Prediction err", "Tolerance", "Leaf Compr.", "Node Compr.");
  }
  else
    f = fopen("Integral.gnu","a");
  
  fprintf(f,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n", TotalTime-RemainingTime, TimeStep, KineticEnergy, ((Detonation)? err:MaxAverage.at(1)),
    PredictionError.at(1), Tolerance, LeafNb*1./(1<<(ScaleNb*Dimension)), NodeNb*1./(1<<(ScaleNb*Dimension)));
  
  fclose(f);
}
/*
_______________________________________________________________________________

Backup data
_______________________________________________________________________________

*/

void AdaptiveMesh::backup()
{
  FILE *f;
  int Nl = N(ScaleNb);
  int l=ScaleNb;
  int i,j,k,n;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  int imin = 0;
  int imax = ei*(Nl-1);
  int jmin = 0;
  int jmax = ej*(Nl-1);
  int kmin = 0;
  int kmax = ek*(Nl-1);
  
  f = fopen("lara.bak","w");
  
  fprintf(f,"%15.8e\n", TimeStep);
  fprintf(f,"%15.8e\n", TotalTime-RemainingTime);
  fprintf(f,"%15i\n", IterationNo);
  fprintf(f,"%15i\n", WriteNo);
 
  for (k=kmin; k<=kmax; k++)
  for (j=jmin; j<=jmax; j++)
  for (i=imin; i<=imax; i++)
  {
    for (n=1;n<=QuantityNb;n++)
      fprintf(f,"%15.8e  ",cell(l,i,j,k)->average(n,1));
    
    fprintf(f,"\n");
  }
  fclose(f);
}

/*
_______________________________________________________________________________

Restore data
_______________________________________________________________________________

*/

void AdaptiveMesh::restore()
{
  FILE *f;
  int Nl = N(ScaleNb);
  int l=ScaleNb;
  int i,j,k,n;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  int imin = 0;
  int imax = ei*(Nl-1);
  int jmin = 0;
  int jmax = ej*(Nl-1);
  int kmin = 0;
  int kmax = ek*(Nl-1);
  float elapsedTime, timeStep;
  float x;
  int rep=0;
  
  f = fopen("lara.bak","r");
  
  if (!f) return;
  
  rep=fscanf(f,"%e", &timeStep);
  rep=fscanf(f,"%e", &elapsedTime);
  rep=fscanf(f,"%i", &IterationNo); 
  rep=fscanf(f,"%i", &WriteNo);
 
  TimeStep = timeStep;
  RemainingTime = TotalTime - elapsedTime;
  StartTime = elapsedTime;
  
  for (k=kmin; k<=kmax; k++)
  for (j=jmin; j<=jmax; j++)
  for (i=imin; i<=imax; i++)
  {
    for (n=1;n<=QuantityNb;n++)
    {
      rep=fscanf(f,"%e  ",&x);
      cell(l,i,j,k)->setAverage(n,1,x);
    }
  }
  
  fclose(f);
}
/*
_______________________________________________________________________________

Penalize
_______________________________________________________________________________

*/

void AdaptiveMesh::penalize()
{ 
  int l, i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;
  
  // --- Penalize momentum for n->n+1/2 ---

  for (l=ScaleMin; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=0; i <= ei*(N(l)-1); i++)
  for (j=0; j <= ej*(N(l)-1); j++)
  for (k=0; k <= ek*(N(l)-1); k++)
  {    
    if (cell(l,i,j,k)->status()==Leaf)
    {
      if (cell(l,i,j,k)->isInSolid()) 
      {
        cell(l,i,j,k)->penalizeMomentum();
	cell(l,i,j,k)->penalizeEnergy();
      }
    }
  }
  
  // --- Penalize density for n->n+1 ---
  
  for (l=ScaleMin; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=0; i <= ei*(N(l)-1); i++)
  for (j=0; j <= ej*(N(l)-1); j++)
  for (k=0; k <= ek*(N(l)-1); k++)
  {    
    if (cell(l,i,j,k)->status()==Leaf)
    {
      if (cell(l,i,j,k)->isInSolid()) 
      {
	cell(l,i,j,k)->setDensity(cell(l,i,j,k)->density()+TimeStep/(4.*cell(l,i,j,k)->size(1))*(cell(l,i-1,j,k)->momentum(1)-cell(l,i+1,j,k)->momentum(1)));

	  if (Dimension > 1)
	    cell(l,i,j,k)->setDensity(cell(l,i,j,k)->density()+TimeStep/(4.*cell(l,i,j,k)->size(2))*(cell(l,i,j-1,k)->momentum(2)-cell(l,i,j+1,k)->momentum(2)));	  
	
	  if (Dimension > 2)
	    cell(l,i,j,k)->setDensity(cell(l,i,j,k)->density()+TimeStep/(4.*cell(l,i,j,k)->size(3))*(cell(l,i,j,k-1)->momentum(3)-cell(l,i,j,k+1)->momentum(3)));	  
/*		
		cell(l,i,j,k)->penalizeDensity(cell(l,i-3,j,k), cell(l,i-2,j,k),cell(l,i-1,j,k), cell(l,i+1,j,k),cell(l,i+2,j,k), cell(l,i+3,j,k),1);
		if (Dimension > 1) 
		  cell(l,i,j,k)->penalizeDensity(cell(l,i,j-3,k), cell(l,i,j-2,k),cell(l,i,j-1,k), cell(l,i,j+1,k),cell(l,i,j+2,k), cell(l,i,j+3,k),2);
		if (Dimension > 2) 
		  cell(l,i,j,k)->penalizeDensity(cell(l,i,j,k-3), cell(l,i,j,k-2),cell(l,i,j,k-1), cell(l,i,j,k+1),cell(l,i,j,k+2), cell(l,i,j,k+3),3);	
*/
      }
    }
  }
  
  // --- Penalize momentum for n+1/2->n+1 ---

  for (l=ScaleMin; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=0; i <= ei*(N(l)-1); i++)
  for (j=0; j <= ej*(N(l)-1); j++)
  for (k=0; k <= ek*(N(l)-1); k++)
  {    
    if (cell(l,i,j,k)->status()==Leaf)
    {
      if (cell(l,i,j,k)->isInSolid()) 
      {
	cell(l,i,j,k)->penalizeMomentum();
 	cell(l,i,j,k)->penalizeEnergy();
      }
   }
  }
  
/*
  
  // --- Penalize energy for n->n+1 ---

  for (l=ScaleMin; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=0; i <= ei*(N(l)-1); i++)
  for (j=0; j <= ej*(N(l)-1); j++)
  for (k=0; k <= ek*(N(l)-1); k++)
  {    
    if (cell(l,i,j,k)->status()==Leaf)
    {
      if (cell(l,i,j,k)->isInSolidAtMid()) 
	cell(l,i,j,k)->penalizeEnergy();
    }
  }
*/
}

/*
_______________________________________________________________________________

Interpolate in cell l,i,j,k when the cell becomes fluid
_______________________________________________________________________________

*/

void AdaptiveMesh::interpolateFluid(const int l, const int i, const int j, const int k)
{
  int n;
  Vector result(QuantityNb);
  
  n = 0;
  
  if (!cell(l,i-1,j,k)->wasInSolid())
  {
    n++;
    result += cell(l,i-1,j,k)->average();
  }
  
  if (!cell(l,i+1,j,k)->wasInSolid())
  {
    n++;
    result += cell(l,i+1,j,k)->average();
  }
  
  if (Dimension > 1 && !cell(l,i,j-1,k)->wasInSolid())
  {
    n++;
    result += cell(l,i,j-1,k)->average();
  }
     
  if (Dimension > 1 && !cell(l,i,j+1,k)->wasInSolid())
  {
    n++;
    result += cell(l,i,j+1,k)->average();
  }
     
   if (Dimension > 2 && !cell(l,i,j,k-1)->wasInSolid())
  {
    n++;
    result += cell(l,i,j,k-1)->average();
  }
     
  if (Dimension > 1 && !cell(l,i,j,k+1)->wasInSolid())
  {
    n++;
    result += cell(l,i,j,k+1)->average();
  }
     
  result *= (1./n);
  
  cell(l,i,j,k)->setAverage(result);
}


/*
_______________________________________________________________________________

Interpolate in cell l,i,j,k when the cell becomes solid
_______________________________________________________________________________

*/

void AdaptiveMesh::interpolateSolid(const int l, const int i, const int j, const int k)
{
  int n;
  Vector result(QuantityNb);
  
  n = 0;
  
  if (cell(l,i-1,j,k)->wasInSolid())
  {
    n++;
    result += cell(l,i-1,j,k)->average();
  }
  
  if (cell(l,i+1,j,k)->wasInSolid())
  {
    n++;
    result += cell(l,i+1,j,k)->average();
  }
  
  if (Dimension > 1 && cell(l,i,j-1,k)->wasInSolid())
  {
    n++;
    result += cell(l,i,j-1,k)->average();
  }
     
  if (Dimension > 1 && cell(l,i,j+1,k)->wasInSolid())
  {
    n++;
    result += cell(l,i,j+1,k)->average();
  }
     
   if (Dimension > 2 && cell(l,i,j,k-1)->wasInSolid())
  {
    n++;
    result += cell(l,i,j,k-1)->average();
  }
     
  if (Dimension > 1 && cell(l,i,j,k+1)->wasInSolid())
  {
    n++;
    result += cell(l,i,j,k+1)->average();
  }
     
  result *= (1./n);
  
  cell(l,i,j,k)->setAverage(result);
}

/*
_______________________________________________________________________________

Render date
_______________________________________________________________________________

*/
void AdaptiveMesh::render()
{
  int l, i, j, k;
  int ei = 1;
  int ej = (Dimension > 1) ? 1:0;
  int ek = (Dimension > 2) ? 1:0;

  if (PenalVelocity == 0.) return;
  
  for (l=ScaleMin; l<=ScaleNb; l++)
#pragma omp parallel for private(i,j,k) shared(l,ei,ej,ek)
  for (i=0; i <= ei*(N(l)-1); i++)
  for (j=0; j <= ej*(N(l)-1); j++)
  for (k=0; k <= ek*(N(l)-1); k++)
  {    
    if (cell(l,i,j,k)->status()==Leaf)
    {
      if(cell(l,i,j,k)->wasInSolid() && !cell(l,i,j,k)->isInSolid())
	interpolateFluid(l,i,j,k);
      
      if(!cell(l,i,j,k)->wasInSolid() && cell(l,i,j,k)->isInSolid())
	interpolateSolid(l,i,j,k);
      
      
      
    }
  }
}

/*
_______________________________________________________________________________

Init mesh
_______________________________________________________________________________

*/

void AdaptiveMesh::init()
{
  cout << "lara: execution\n";
  
  initAverage();
  if (Restart) restore();
  computeIntegralValues();
  updateAll();
  setBoundaryConditions();
  adapt();
  updateAll();
  setBoundaryConditions();  
  countLeavesAndCells();
  
  if (!Restart)
  {
    writeTree();
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

void AdaptiveMesh::end()
{
  writeTree();
  writeAverage();
  if (WriteStrip) writeStrip();
  computeIntegralValues();
  writeIntegral();
  showPerformance();
}

/*
_______________________________________________________________________________

Compute complete iteration
_______________________________________________________________________________

*/

void AdaptiveMesh::iterate()
{
  // Init iteration
  
  computeTimeStep();
  
  // For detonation compute source
  
  if (Detonation) 
  {  
    storeAverage();
    computeSource();
    updateAll();
    setBoundaryConditions();
  }
   
  if (UsePenalization) 
  {
    render();
    penalize();
  }
   
  // Multi-stage method

  storeAverage();
  
  for (StepNo = 1; StepNo <= StepNb; StepNo++)
  {
    if (Equation == NavierStokes && !IsInviscid) 
      computeGradient();
    
    computeDivergence();
    computeStep();
    updateAll();
    setBoundaryConditions();
  }

  // For detonation compute source
  
  if (Detonation) 
  {  
    storeAverage();
    computeSource();
    updateAll();
    setBoundaryConditions();
  }

  if (UsePenalization) penalize();

  // Compute integral values
  
  checkStability();
  computeIntegralValues();

  // Adapt tree structure
  
  adapt();
  updateAll();
  setBoundaryConditions();
  countLeavesAndCells();
  
  // Write cell-average values into file if required
  
  if (timeToWrite())
  {
     writeTree();
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
