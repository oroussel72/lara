#include <iostream>
#include <stdio.h>
#include <stdlib.h>
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

/*
_______________________________________________________________________________

Constructor
_______________________________________________________________________________

*/

Cell::Cell() : 
  _center(Dimension),
  _size(Dimension),
  _average(QuantityNb,StorageNb),
  _gradient(Dimension,Dimension)
{
  ;
}

/*
_______________________________________________________________________________

Distructor
_______________________________________________________________________________

*/

Cell::~Cell()
{
}

/*
_______________________________________________________________________________

Set average
_______________________________________________________________________________

*/

void Cell::setAverage(const int storageNo, const Vector& V)
{
  int quantityNo;
  
  for (quantityNo=1;quantityNo <= QuantityNb; quantityNo++)
    setAverage(quantityNo, storageNo, V.at(quantityNo));  
}


/*
_______________________________________________________________________________

Reset average
_______________________________________________________________________________

*/

void Cell::resetAverage(const int storageNo)
{
  int quantityNo;
  
  for (quantityNo=1;quantityNo <= QuantityNb; quantityNo++)
    setAverage(quantityNo, storageNo, 0.);  
}

/*
_______________________________________________________________________________

Velocity vector
_______________________________________________________________________________

*/

Vector Cell::velocity(const int storageNo) const
{
	Vector V(Dimension);
	int axisNo;

	for (axisNo = 1; axisNo <= Dimension; axisNo++)
		V.set(axisNo, velocity(axisNo,storageNo));

	return V;
}


/*
_______________________________________________________________________________

Velocity norm
_______________________________________________________________________________

*/

real Cell::velocityNorm(const int storageNo) const
{
  switch (Dimension)
  {
    case 1:
    default:
      return ABS(velocity(1,storageNo));
      break;
      
    case 2:
      return sqrt(SQR(velocity(1,storageNo))+SQR(velocity(2,storageNo)));
      break;
      
    case 3:
      return sqrt(SQR(velocity(1,storageNo))+SQR(velocity(2,storageNo))+SQR(velocity(3,storageNo)));
      break;
  };
}

/*
_______________________________________________________________________________

Momentum norm
_______________________________________________________________________________

*/

real Cell::momentumNorm(const int storageNo) const
{
  switch (Dimension)
  {
    case 1:
    default:
      return ABS(momentum(1,storageNo));
      break;
      
    case 2:
      return sqrt(SQR(momentum(1,storageNo))+SQR(momentum(2,storageNo)));
      break;
      
    case 3:
      return sqrt(SQR(momentum(1,storageNo))+SQR(momentum(2,storageNo))+SQR(momentum(3,storageNo)));
      break;
  };
}
/*
_______________________________________________________________________________

Maximal eigenvalue for the time step computation
_______________________________________________________________________________

*/

real Cell::maxEigenvalue(const int storageNo) const
{
  switch(Equation)
  {
    case ConvectionDiffusion:
    default:      
      return ABS(Celerity);
      break;
      
    case Burgers:
      return velocityNorm(storageNo);
      break;
      
    case NavierStokes:
      return velocityNorm(storageNo)+speedOfSound(storageNo);
      break;
  };
}

/*
_______________________________________________________________________________

Pressure
_______________________________________________________________________________

*/

real Cell::pressure(const int storageNo) const
{
	real rho  = density(storageNo);
	real rhoE = energy(storageNo);
	Vector V  = velocity(storageNo);
	real rhoY = (SpeciesNb > 1) ? average(Dimension+3,storageNo):0.;

	return EOS_P(rho, V, rhoE, rhoY);	
}

/*
_______________________________________________________________________________

Temperature
_______________________________________________________________________________

*/

real Cell::temperature(const int storageNo) const
{
	real rho = density(storageNo);
	real p = pressure(storageNo);
	real M = MolarMass;

	return M*p/(rho*R);
}

/*
_______________________________________________________________________________

Return true if a nan is detected
_______________________________________________________________________________

*/

bool Cell::isNaN() const
{
  int i;
  
  for (i=1; i<=QuantityNb; i++)
    if (isnan(average(i,1))) return true;
    
  return false;
}


/*
_______________________________________________________________________________

Compute inviscid Euler flux
_______________________________________________________________________________

*/

Vector Cell::eulerFlux(const Vector& average, const int axisNo) const
{
	int i;
	real rho, p, rhoE, u, rhoY;
	Vector result(QuantityNb);
	Vector V(Dimension);
	
	rho=p=rhoE=u=rhoY=0.;
	
	switch (Equation)
	{
	  case ConvectionDiffusion:
	  default:
	    u = average.at(1);
	    result.set(1, Celerity*u);
	    break;
	    
	  case Burgers:
	    for (i = 1; i <= Dimension; i++)
	      result.set(i, average.at(i)*average.at(axisNo));
	    break;

	  case NavierStokes:
	    rho = average.at(1);
	    
	    for (i=1; i<= Dimension; i++)
	      V.set(i,average.at(i+1)/rho);
	    
	    rhoE = average.at(Dimension+2);
	    
	    if (SpeciesNb > 1)
	      rhoY = average.at(Dimension+3);
	    
	    p = EOS_P(rho, V, rhoE, rhoY);
	    
	    result.set(1, rho*V.at(axisNo));
	
	    for (i = 1; i <= Dimension; i++)
	      result.set(i+1, rho*V.at(i)*V.at(axisNo) + ((axisNo==i) ? p : 0.) ); 
	    
	    result.set(Dimension+2, (rhoE + p)*V.at(axisNo));

	    if (SpeciesNb > 1)
	    {
	      for (i = 1; i<= SpeciesNb-1; i++)
		result.set(Dimension+2+i, average.at(Dimension+2+i)*V.at(axisNo));
	    }
	    
	    if (UsePenalization)
	      result.set(QuantityNb, 0.);
	    
	    break;
	};
	
	return result;
}


/*
_______________________________________________________________________________

Copy stage sourceStagoNo into targetStageNo
_______________________________________________________________________________

*/


void Cell::copyStorage(const int sourceStorageNo, const int targetStorageNo)
{
	int quantityNo;

	for (quantityNo = 1; quantityNo <= QuantityNb; quantityNo++)
		setAverage(quantityNo,targetStorageNo, average(quantityNo,sourceStorageNo));
}

/*
_______________________________________________________________________________

Compute velocity gradient
_______________________________________________________________________________

*/

void Cell::computeGradient(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const int axisNo)
{
  // Only for Navier-Stokes

  if (Equation != NavierStokes) return;

  int i;
  real result;
  real dx = size(axisNo);
  real a = 1./(12.*dx);
  
  // Compute gradient

  for (i = 1; i <= Dimension; i++)
  {
    result = a*(-Cip2->v(i)+8.*(Cip1->v(i)-Cim1->v(i))+Cim2->v(i));
    setGradient(i, axisNo, result);
  }
}
/*
_______________________________________________________________________________

Inviscid flux
_______________________________________________________________________________

*/

Vector Cell::inviscidFlux(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo) const
{
	switch(Scheme)
	{
		case McCormack:
		default:
			return schemeMcCormack(Cim1, Cip1, Cip2, axisNo);		
			break;
		
		case Roe:
			return schemeRoe(Cim2, Cim1, Cip1, Cip2, Cip3, axisNo);
			break;
			
		case AUSM:
 			return schemeAUSM(Cim2, Cim1, Cip1, Cip2, Cip3, axisNo);
			break;
	};

}

/*
_______________________________________________________________________________

2-4 McCormack scheme in (i + 1/2) [i p 1/2]
_______________________________________________________________________________

*/

Vector Cell::schemeMcCormack(const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const int axisNo) const
{
	int direction; 				// direction (upwind or downwind)
	int loopIterationNo;			// iteration number in the local loop

	Vector average1(QuantityNb);		// Cell-average values
	Vector average2(QuantityNb);		//
	
	
	Vector Flux1(QuantityNb);		// Fluxes
	Vector Flux2(QuantityNb);		//

	// Compute direction

	loopIterationNo = IterationNo%(1<<Dimension);

	if (axisNo == 1)
		direction = loopIterationNo%2;
	else if (axisNo == 2)
		direction = (loopIterationNo/2)%2;
	else
		direction = (loopIterationNo/2)/2;
	
	direction = (direction + (StepNo-1))%2;

	// Following direction, choose cells

	if (direction == 0)
	{
	  average1 = Cim1->average();
	  average2 = average();
	}
	else
	{
	  average1 = Cip2->average();
	  average2 = Cip1->average();
	}

	// Compute inviscid fluxes

	Flux1 = eulerFlux(average1, axisNo);
	Flux2 = eulerFlux(average2, axisNo);

	// 2-4 Mc Comrmack scheme

	return (7./6.)*Flux2 - (1./6.)*Flux1;
}

/*
_______________________________________________________________________________

Lax-Friedrichs scheme for density (only for penalization)
_______________________________________________________________________________

*/

real Cell::schemeRoeDensity(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo) const
{
  Vector averageL(1);  // Left average
  Vector averageR(1);  // Right average

  Vector Aim2(1);   	// 
  Vector Aim1(1);	//
  Vector Ai(1);		// Cell-average values for density and momentum
  Vector Aip1(1);	//
  Vector Aip2(1);	//
  Vector Aip3(1);	//
  
  real rhoUL, rhoUR, rhoUC;
    
  real result;
    
  int  WENOVariables2 = WENOVariables;	// Store WENOVariables
  
  // Impose conservative WENO variables
  
  WENOVariables = Conservative;
  
  // Store momentum in the direction axisNo

  Aim2.set(1,Cim2->average(axisNo+1,1));
  Aim1.set(1,Cim1->average(axisNo+1,1));
  Ai.  set(1,average(axisNo+1,1));
  Aip1.set(1,Cip1->average(axisNo+1,1));
  Aip2.set(1,Cip2->average(axisNo+1,1));
  Aip3.set(1,Cip3->average(axisNo+1,1));  
  
  
  // Compute WENO reconstruction
  
  averageL = WENOAverage(Aim2, Aim1, Ai, Aip1, Aip2,axisNo);
  averageR = WENOAverage(Aip3, Aip2, Aip1, Ai, Aim1,axisNo);
  
  rhoUL = averageL.at(1);
  rhoUR = averageR.at(1);
  rhoUC = 0.5*(rhoUL+rhoUR);
  
  // Roe scheme (ODF formulation)
  
  result = 0.5*((rhoUL+rhoUR) - SGN(rhoUC)*(rhoUR-rhoUL));
  
  // Turn back to WENOVariables
  
  WENOVariables = WENOVariables2;
  
  return result;
}


/*
_______________________________________________________________________________

Roe scheme
_______________________________________________________________________________

*/

Vector Cell::schemeRoe(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo) const
{
  Vector averageL(QuantityNb);  // Left average
  Vector averageR(QuantityNb);  // Right average
  Vector averageC(QuantityNb);  // Centered average with Roe formula
  Vector result(QuantityNb);
  
  Matrix A(QuantityNb, QuantityNb); 	// Jacobian matrix
  Matrix D(QuantityNb, QuantityNb); 	// Diagonal matrix D with absolute eigenvalues
  Matrix L(QuantityNb, QuantityNb); 	// Left eigenmatrix
  Matrix R(QuantityNb, QuantityNb); 	// Right eigenmatrix
  int i;
      
  switch(Reconstruction)
  {
    case None:
      averageL = average();
      averageR = Cip1->average();
    
    case MUSCL:
      averageL = leftLimitedAverage(Cim1->average(),average(),Cip1->average());
      averageR = rightLimitedAverage(average(),Cip1->average(),Cip2->average());
      break;
      
    case WENO:
      averageL = WENOAverage(Cim2->average(), Cim1->average(), average(), Cip1->average(), Cip2->average(),axisNo);
      averageR = WENOAverage(Cip3->average(), Cip2->average(), Cip1->average(), average(), Cim1->average(),axisNo);
      break;
  }; 
         
  // Compute Roe average state
  
  averageC = roeAverage(averageL, averageR);

  switch(Equation)
  {
    case ConvectionDiffusion:
      A.set(1,1,ABS(Celerity));
      break;
      
    case Burgers:
      for (i=1; i <= Dimension; i++)
	A.set(i,i,ABS(averageC.at(i)));
	break;
    
    case NavierStokes:
      
      // Compute Jacobian matrix
  
      D = absDiagonalMatrix(averageC,axisNo);
      L = leftEigenMatrix(averageC,axisNo);
      R = rightEigenMatrix(averageC,axisNo);
  
      A = R*D*L;
      break;
  };

  // Compute Roe scheme
    
  result = 0.5*(eulerFlux(averageL,axisNo)+eulerFlux(averageR,axisNo)-A*(averageR-averageL));

//  if (UsePenalization && isInSolid()) result.set(1, (1./Porosity)*result.at(1));
  
  return result;
}

/*
_______________________________________________________________________________

AUSM+ scheme
_______________________________________________________________________________

*/
Vector Cell::schemeAUSM(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo) const
{
  Vector averageL(QuantityNb);  // Left average
  Vector averageR(QuantityNb);  // Right average
  Vector result(QuantityNb);
  
  real rhoL, rhoR, rho;		// densities (left, right, center)
  real pL, pR, p;		// pressures
  real rhoEL, rhoER;		// energies
  real rhoHL, rhoHR, rhoH;	// enthalpies
  real cL, cR, c;		// speeds of sound
  Vector VL(Dimension);		//
  Vector VR(Dimension);		// velocities
  Vector V(Dimension);		//
  real ML, MR, M;		// Right, left and central Mach numbers
  real MLp, MRm; 		// AUSM variables
  real PLp, PRm;		// idem
  int i;			// counter
  real M1, M2;			// intermediate computations
  Vector rhoY(SpeciesNb);	//
  Vector rhoYL(SpeciesNb);	// Species
  Vector rhoYR(SpeciesNb);	//
  
  real a = 3./16.;		// AUSM parameters
  real b = 0.125;		//
  
  if (Equation != NavierStokes)
  {  
    printf("lara: in function Vector Cell::schemeAUSM(const Cell*, const Cell*, const Cell*, const Cell*, const Cell*, const int) const\n");
    printf("lara: only for Equation = NavierStokes\n");
    exit(-1);
  }
    
  switch(Reconstruction)
  {
    case None:
      averageL = average();
      averageR = Cip1->average();
    
    case MUSCL:
      averageL = leftLimitedAverage(Cim1->average(),average(),Cip1->average());
      averageR = rightLimitedAverage(average(),Cip1->average(),Cip2->average());
      break;
      
    case WENO:
      averageL = WENOAverage(Cim2->average(), Cim1->average(), average(), Cip1->average(), Cip2->average(),axisNo);
      averageR = WENOAverage(Cip3->average(), Cip2->average(), Cip1->average(), average(), Cim1->average(),axisNo);
      break;
  }; 
  
  // --- Compute primitive variables ---
  
  rhoL = averageL.at(1);
  rhoR = averageR.at(1);
  
  for (i=1; i <= Dimension; i++)
  {
    VL.set(i,averageL.at(i+1)/rhoL);
    VR.set(i,averageR.at(i+1)/rhoR);
  }
  
  rhoEL = averageL.at(Dimension+2);
  rhoER = averageR.at(Dimension+2);
  
  if (SpeciesNb > 1)
  {
  
    for (i=1; i <= SpeciesNb-1; i++)
    {
      rhoYL.set(i,averageL.at(Dimension+2+i));
      rhoYR.set(i,averageR.at(Dimension+2+i));      
    }
  }
  
  pL = EOS_P(rhoL, VL, rhoEL, rhoYL.at(1));
  pR = EOS_P(rhoR, VR, rhoER, rhoYR.at(1));
  
  rhoHL = rhoEL + pL;
  rhoHR = rhoER + pR;
  
  cL = sqrt(Gamma*pL/rhoL);
  cR = sqrt(Gamma*pR/rhoR);
  
  c = sqrt(cL*cR);
    
  ML = VL.at(axisNo)/c;
  MR = VR.at(axisNo)/c;
  
  // --- Compute AUSM variables ---
    
  // Compute MLp [ M+(i) ],MRm [ M-(i+1) ], PLp [ P+(i) ] PRm [ P-(i+1) ]
  
  if (ABS(ML) >= 1.)
  {
    MLp = 0.5*(ML+ABS(ML));
    PLp = 0.5*(1.+SGN(ML));
  }
  else
  {
    M1 = SQR(ML+1.);
    M2 = SQR(ML*ML-1.);
    
    MLp = 0.25*M1 + b*M2;
    PLp = 0.25*M1*(2.-ML) + a*ML*M2;
  }
  
  if (ABS(MR) >= 1.)
  {
    MRm = 0.5*(MR-ABS(MR));
    PRm = 0.5*(1.-SGN(ML));
  }
  else
  {
    M1 = SQR(MR-1.);
    M2 = SQR(MR*MR-1.);
    
    MRm = -0.25*M1 - b*M2;
    PRm = 0.25*M1*(2.+MR) - a*MR*M2;
  }
  
  // Compute p and M
    
  M = MLp + MRm;
  p = PLp*pL + PRm*pR;
    
  // Decenter flux
  
  if (M > 0.)
  {
    rho = rhoL;
    V = VL;
    rhoH = rhoHL;
    if (SpeciesNb > 1) rhoY = rhoYL;
    
  }
  else if (M < 0.)
  {
    rho = rhoR;
    V = VR;
    rhoH = rhoHR;
    if (SpeciesNb > 1) rhoY = rhoYR;
  }
  else
  {
    rho = 0.5*(rhoL+rhoR);
    V = 0.5*(VL+VR);
    rhoH = 0.5*(rhoHL+rhoHR);
    if (SpeciesNb > 1) rhoY = 0.5*(rhoYL+rhoYR);
  }
  
  result.set(1,rho);
  for(i=1; i<= Dimension; i++)
    result.set(i+1,rho*V.at(i));
  result.set(Dimension+2,rhoH);
  
  if (SpeciesNb > 1)
  {
    for (i=1; i<= SpeciesNb-1; i++)
      result.set(Dimension+2+i,rhoY.at(i));
  }
    
  result *= M*c;
  
  result.set(1+axisNo,result.at(1+axisNo)+p);
    
  return result;
}

/*
_______________________________________________________________________________

Viscous flux
_______________________________________________________________________________

*/
Vector Cell::viscousFlux(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo) const
{
  Vector result(QuantityNb);
  Vector tau(Dimension);
  Vector V(Dimension);
  real dx, mu, temperature;
  Matrix gradV(Dimension,Dimension);
  real gradT, divV;
  Vector gradY(SpeciesNb);
  int i,j;
  real a = 1./12.;
  real b;
  real twoThird = 2./3.;
  
  dx = 0.5*(size(axisNo)+Cip1->size(axisNo));
  b = 1./(12.*dx);
  
  if (Equation != NavierStokes)
    result = -b*Diffusivity*(Cim1->average()+15.*(-average()+Cip1->average())-Cip2->average());
  else
  {       
    // --- Compute interface values in i+1/2 ---
    
    // Compute gradV
    
    for (i=1; i<=Dimension; i++)
    for (j=1; j<=Dimension; j++)
    {
      if (j != axisNo)
	gradV.set(i,j, a*(-Cim1->gradient(i,j)+7.*(gradient(i,j)+Cip1->gradient(i,j))-Cip2->gradient(i,j)));
    }
    
    for (i=1; i<=Dimension; i++)
      gradV.set(i, axisNo, b*(Cim1->v(i)+15.*(-v(i)+Cip1->v(i))-Cip2->v(i)));
    
    // Compute divV
    
    divV = 0.;
    for (i=1; i<=Dimension; i++)
      divV += gradV.at(i,i);
    
    // Compute V
      
    for (i=1; i<=Dimension; i++)
      V.set(i, a*(-Cim1->v(i)+7.*(v(i)+Cip1->v(i))-Cip2->v(i)));

    //  Compute temperature
      
    temperature = a*(-Cim1->T()+7*(T()+Cip1->T())-Cip2->T());
    
    // Compute gradT
      
    gradT = b*(Cim1->T()+15.*(-T()+Cip1->T())-Cip2->T());   

    // Compute gradY
    
    if (SpeciesNb > 1)
    {
      for (i=1; i <= SpeciesNb-1; i++)
	gradY.set(i,b*(Cim1->Y(i)+15.*(-Y(i)+Cip1->Y(i))-Cip2->Y(i)));
    }
    
    // Compute viscosity
    
    mu = sutherland(temperature);   
    
    // Compute tau
   
    for (i=1; i<=Dimension; i++)
      tau.set(i,mu*(gradV.at(axisNo,i) + gradV.at(i,axisNo) - ( (axisNo==i)? twoThird*divV : 0. )));
    
    // Compute flux
      
    result.set(1,0.);
    
    for (i=1; i<= Dimension; i++)
      result.set(i+1, -tau.at(i));
    
    result.set(Dimension+2, -tau*V - mu*Conductivity*gradT - ((Detonation)? Q0*MassDiffusivity.at(1)*gradY.at(1):0.));
    
    if (SpeciesNb > 1)
    {
      for (i=1; i <= SpeciesNb-1; i++)
	result.set(Dimension+2+i, -MassDiffusivity.at(i)*gradY.at(i));
    }
  }
      
  return result;
}

/*
_______________________________________________________________________________

Limiter
_______________________________________________________________________________

*/

real Cell::limiter(const real r) const
{
  switch(Limiter)
  {
    case None:
    default:
      return 0.;
    
    case VanAlbada:
      return MAX(0.,(SQR(r)+r)/(SQR(r)+1.));
      
    case VanLeer:
      return (r+ABS(r))/(1.+ABS(r));
      
    case Koren:
      return MAX(0.,MIN(2.*r, MIN((1.+2.*r)/3.,2.)));
      
    case Ospre:
      return MAX(0.,(3.*(SQR(r)+r))/(2.*(SQR(r)+r+1.)));
      
    case MinMod:
      return MAX(0.,MIN(1.,r));
      
    case SuperBee:
      return MAX(0.,MAX(MIN(1.,2.*r),MIN(2.,r)));
    
  };
}
  

/*
_______________________________________________________________________________

Compute Roe average
_______________________________________________________________________________

*/

Vector Cell::roeAverage(const Vector& average1, const Vector& average2) const
{
  real rho1, rho2, rho;
  real rhoE1, rhoE2, rhoE;
  real rhoH1, rhoH2, rhoH;
  real p1, p2;
  Vector rhoY(SpeciesNb);
  Vector rhoY1(SpeciesNb);
  Vector rhoY2(SpeciesNb);
  Vector V1(Dimension);
  Vector V2(Dimension);
  Vector V(Dimension);
  Vector result(QuantityNb);
  real a, a1, a2;
  int i;
  
  switch(Equation)
  {
    // For Convection-diffusion and Burgers, apply a standard average
    
    case ConvectionDiffusion:
    case Burgers:
    default:
      result = 0.5*(average1+average2);
      break;
      
    case NavierStokes:
      
      // Compute primitive variables
      
      rho1 = average1.at(1);
      rho2 = average2.at(1);

      for (i = 1; i <= Dimension; i++)
      {
	V1.set(i,average1.at(i+1)/rho1);
	V2.set(i,average2.at(i+1)/rho2);
      }
      
      rhoE1 = average1.at(Dimension+2);
      rhoE2 = average2.at(Dimension+2);
      
      if (SpeciesNb > 1)
      {
	for (i=1; i <= SpeciesNb-1; i++)
	{
	  rhoY1.set(i, average1.at(Dimension+2+i));
	  rhoY2.set(i, average2.at(Dimension+2+i));
	}
      }
      
      p1 = EOS_P(rho1, V1, rhoE1, rhoY1.at(1));
      p2 = EOS_P(rho2, V2, rhoE2, rhoY2.at(1));
      
      rhoH1 = rhoE1 + p1;
      rhoH2 = rhoE2 + p2;
      
      // Compute Roe average values
      
      a = sqrt(rho2/rho1);
      a1 = (1./(1. + a));
      a2 = 1.-a1;
      
      rho = a*rho1;
      V = a1*V1 + a2*V2;
      rhoH = a1*rhoH1 + a2*rhoH2;
      rhoY = a1*rhoY1 + a2*rhoY2;
      
      // Compute cell-average vector;
      
      rhoE = 1./Gamma*(rhoH + (Gamma-1.)*(0.5*rho*V*V + ((Detonation)? rhoY.at(1)*Q0:0.)));
      
      result.set(1,rho);
      
      for (i=1;i<=Dimension;i++)
	result.set(i+1,rho*V.at(i));
       
      result.set(Dimension+2,rhoE);
      
      if (SpeciesNb > 1)
      {
	for (i=1; i <= SpeciesNb-1; i++)
	  result.set(Dimension+2+i,rhoY.at(i));
      }
      
      break;
      
  };
  
  return result; 
}

/*
_______________________________________________________________________________

Left eigenmatrix
_______________________________________________________________________________

*/

Matrix Cell::leftEigenMatrix(const Vector& average, const int axisNo) const
{

  // --- Local variables ---

  Vector V(Dimension);
  real rho, rhoE, p, c;
  real V2; 			// square of the velocity norm
  real h; 			// internal enthalpy
  int i,j;

  real invH, invC;
  
  Matrix result(QuantityNb,QuantityNb);

  // --- Compute local variables ---
  
  rho = average.at(1);
    
  for (i=1; i<= Dimension; i++)
    V.set(i,average.at(i+1)/rho);
  
  rhoE = average.at(Dimension+2);
  V2 = V*V;
  p = (Gamma-1.)*(rhoE - 0.5*rho*V2);
  c = sqrt(Gamma*p/rho);
  h = c*c/(Gamma - 1.);
  invH = 1./h;
  invC = 1./c;
  
  
  // --- 1st column ---
  
  for (i=1; i<=Dimension;i++)
  {
    if (i==axisNo)
      result.set(i,1,1. - 0.5*invH*V2);
    else
      result.set(i,1,-V.at(i));
  }
  
  result.set(Dimension+1, 1, 0.5*( 0.5*invH*V2 - V.at(axisNo)*invC));
  result.set(Dimension+2, 1, 0.5*( 0.5*invH*V2 + V.at(axisNo)*invC));
  
  
  // --- 2nd (3rd-4th) columns ---
  
  for (i=1; i <= Dimension; i++)
  for (j=1; j <= Dimension; j++)
  {
    if (i==axisNo)
      result.set(i,j+1,V.at(j)*invH);
    else if (i==j)
      result.set(i, j+1, 1.);
    else
      result.set(i, j+1, 0.);
  }
  
  for (j=1; j <= Dimension; j++)
  {
    if (j==axisNo)
    {
      result.set(Dimension+1, j+1, 0.5*( invC-V.at(axisNo)*invH));
      result.set(Dimension+2, j+1, 0.5*(-invC-V.at(axisNo)*invH));
    }
    else
    {
      result.set(Dimension+1, j+1, -0.5*V.at(j)*invH);
      result.set(Dimension+2, j+1, -0.5*V.at(j)*invH);
    }
  }
  
  // Last columns
  
  for (i=1; i <= Dimension; i++)
  {
    if (i==axisNo)
      result.set(i,Dimension+2,-invH);
    else
      result.set(i,Dimension+2, 0.);
  }
  
  result.set(Dimension+1, Dimension+2, 0.5*invH);
  result.set(Dimension+2, Dimension+2, 0.5*invH);
    
  return result;
}

/*
_______________________________________________________________________________

Right eigenmatrix
_______________________________________________________________________________

*/

Matrix Cell::rightEigenMatrix(const Vector& average, const int axisNo) const
{

  // --- Local variables ---

  Vector V(Dimension);
  real rho, rhoE, p, c, V2;
  real h; 			// total enthalpy
  int i,j;
  
  Matrix result(QuantityNb,QuantityNb);

  // --- Compute local variables ---
  
  rho = average.at(1);
  
  for (i=1; i<= Dimension; i++)
    V.set(i,average.at(i+1)/rho);
  
  rhoE = average.at(Dimension+2);
  
  V2 = V*V;
  p = (Gamma-1.)*(rhoE - 0.5*rho*V2);
  c = sqrt(Gamma*p/rho);
  h = (rhoE +p)/rho;

  // --- First line ---
  
  for (j=1; j <= Dimension; j++)
    result.set(1,j,(j==axisNo)? 1.:0.);
  
  result.set(1, Dimension+1, 1.);
  result.set(1, Dimension+2, 1.);
  
  // --- 2nd (3rd-4th) lines ---
  
  for (i=1; i <= Dimension; i++)
  for (j=1; j <= Dimension; j++)
  {
    if (j==axisNo)
      result.set(i+1, j, V.at(i));
    else if (i==j)
      result.set(i+1, j, 1.);
    else
      result.set(i+1, j, 0.);
  }

  for (i=1; i <= Dimension; i++)
  {
    result.set(i+1, Dimension+1, V.at(i) + ((i==axisNo) ? c:0.));
    result.set(i+1, Dimension+2, V.at(i) - ((i==axisNo) ? c:0.));
  }
  
  // --- Last line ---
  
  for (j=1;j <= Dimension; j++)
  {
    if (j==axisNo)
      result.set(Dimension+2, j, 0.5*V2);
    else
      result.set(Dimension+2, j, V.at(j));
  }
  
  result.set(Dimension+2, Dimension+1, h + V.at(axisNo)*c);
  result.set(Dimension+2, Dimension+2, h - V.at(axisNo)*c);
  
  return result;
}

/*
_______________________________________________________________________________

Diagonal matrix D containing the absolute eigenvalues
_______________________________________________________________________________

*/

Matrix Cell::absDiagonalMatrix(const Vector& average, const int axisNo) const
{
  // --- Local variables ---

  Vector V(Dimension);
  real rho, rhoE, p, c;
  int i;
  
  Matrix result(QuantityNb,QuantityNb);

  // --- Compute local variables ---
  
  rho = average.at(1);
  
  for (i=1; i<= Dimension; i++)
    V.set(i,average.at(i+1)/rho);
  
  rhoE = average.at(Dimension+2);
  
  p = (Gamma-1.)*(rhoE - 0.5*rho*V*V);
  c = sqrt(Gamma*p/rho);
  
  for (i = 1; i <= Dimension; i++)
    result.set(i, i, ABS(V.at(axisNo)));
  
  result.set(Dimension+1, Dimension+1, ABS(V.at(axisNo)+c));
  result.set(Dimension+2, Dimension+2, ABS(V.at(axisNo)-c));
  
  return result;
}

/*
_______________________________________________________________________________

Left average value with limiter
_______________________________________________________________________________

*/

Vector Cell::leftLimitedAverage(const Vector& averageIm1, const Vector& averageI, const Vector& averageIp1) const
{
  int n;
  real psi;
  real r;
  real eps = 1.E-99;
  Vector result(QuantityNb);
  
  // Compute psi
  
  switch (Centering)
  {
    case Decentered:
    default:
      psi = -1.;
      break;
      
    case Centered:
      psi = 1.;
      break;
      
    case ThirdOrder:
      psi = 1./3.;
      break;
  };
  
  
  for (n=1; n<=QuantityNb; n++)
  {
      // --- Compute average value in i+1/2 on the left side ---
	
      r = (averageIp1.at(n)-averageI.at(n)+eps)/(averageI.at(n)-averageIm1.at(n)+eps);
      result.set(n,averageI.at(n) + 0.25*(1.-psi)*limiter(r   )*(averageI.at(n)-averageIm1.at(n))
				   + 0.25*(1.+psi)*limiter(1./r)*(averageIp1.at(n)-averageI.at(n)));
  }
  
  return result;
}

/*
_______________________________________________________________________________

Right average value with limiter
_______________________________________________________________________________

*/

Vector Cell::rightLimitedAverage(const Vector& averageI, const Vector& averageIp1, const Vector& averageIp2) const
{
  int n;
  real psi;
  real r;
  real eps = 1.E-99;
  Vector result(QuantityNb);
  
  // Compute psi
  
  switch (Centering)
  {
    case Decentered:
    default:
      psi = -1.;
      break;
      
    case Centered:
      psi = 1.;
      break;
      
    case ThirdOrder:
      psi = 1./3.;
      break;
  };
  
  
  for (n=1; n<=QuantityNb; n++)
  {
      // --- Compute average value in i+1/2 on the left side ---
	
      r = (averageIp2.at(n)-averageIp1.at(n)+eps)/(averageIp1.at(n)-averageI.at(n)+eps);
      result.set(n,averageIp1.at(n) - 0.25*(1.+psi)*limiter(r   )*(averageIp1.at(n)-averageI.at(n))
				     - 0.25*(1.-psi)*limiter(1./r)*(averageIp2.at(n)-averageIp1.at(n)));
  }
  
  return result;
}

/*
_______________________________________________________________________________

WENO reconstruction vector
_______________________________________________________________________________

*/

Vector Cell::WENOAverage(const Vector& averageIm2, const Vector& averageIm1, const Vector& averageI, const Vector& averageIp1, const Vector& averageIp2, const int axisNo) const
{
  real g1 = 0.1;		//
  real g2 = 0.6;		// optimal interpolation coefficients
  real g3 = 0.3;		//
  real eps = 1.E-06;		// parameters to avoid zero denominator
  real IS1, IS2, IS3;		// ponderated residuals
  real w1, w2, w3;		// weight coefficients
  real a1, a2, a3, a;		
  
  int i;
  real b = 1./12.;
  real c = 2*b;
  int Nl = averageIm2.size();
    
  Vector U1, U2, U3;		// 3rd order interpolations of the cell-average values
  Vector UIm2, UIm1, UI, UIp1, UIp2; // Average vector after the choice of the variables: characteristic, conservative, primitive
  Vector averageC(QuantityNb);	// Roe average vector
  Matrix L;			// Left eigenmatrix

  Vector result (Nl);
  
  // Compute UIm2, UIm1, UI, UIp1, UIp2
  
  if (Equation != NavierStokes) 
    WENOVariables = Conservative;
  
  switch (WENOVariables)
  {
    case Conservative:
    default:
      UIm2 = averageIm2;
      UIm1 = averageIm1;
      UI   = averageI;
      UIp1 = averageIp1;
      UIp2 = averageIp2;
      break;
      
    case Primitive:
      UIm2 = conservativeToPrimitive(averageIm2);
      UIm1 = conservativeToPrimitive(averageIm1);
      UI   = conservativeToPrimitive(averageI  );
      UIp1 = conservativeToPrimitive(averageIp1);
      UIp2 = conservativeToPrimitive(averageIp2);
      break;
   
    case Characteristic:
      averageC = roeAverage(averageI,averageIp1);
      L = leftEigenMatrix(averageC,axisNo);
      
      UIm2 = L*averageIm2;
      UIm1 = L*averageIm1;
      UI   = L*averageI  ;
      UIp1 = L*averageIp1;
      UIp2 = L*averageIp2;
      break;
  };
  
  
  // Compute interpolations

  U1 = c*(2.*UIm2 - 7.*UIm1 + 11.*UI);
  U2 = c*(-UIm1 + 5.*UI + 2.*UIp1);
  U3 = c*(2.*UI + 5.*UIp1 - UIp2);

  // Compute divided differences
  
  for (i=1; i <= Nl; i++)
  {
    IS1 = b*( 13.*SQR(UIm2.at(i)-2.*UIm1.at(i)+UI.at(i)) +
		    3.*SQR(UIm2.at(i)-4.*UIm1.at(i)+3.*UI.at(i)));
    IS2 = b*( 13.*SQR(UIm1.at(i)-2.*UI.at(i)+UIp1.at(i)) +
		    3.*SQR(UIm1.at(i)-4.*UI.at(i)+3.*UIp1.at(i)));
    IS3 = b*( 13.*SQR(UI.at(i)-2.*UIp1.at(i)+UIp2.at(i)) +
		    3.*SQR(UI.at(i)-4.*UIp1.at(i)+3.*UIp2.at(i)));

    // Compute a1, a2, a3, a;
  
    a1 = g1/SQR(eps+IS1);
    a2 = g2/SQR(eps+IS2);
    a3 = g3/SQR(eps+IS3);
    a = 1./(a1+a2+a3);
  
    // Compute weight coefficients
  
    w1 = a1*a;
    w2 = a2*a;
    w3 = a3*a;
  
    // Compute WENO reconstruction
  
    result.set(i, w1*U1.at(i) + w2*U2.at(i) + w3*U3.at(i));
  }
  
  switch (WENOVariables)
  {
    case Conservative:
    default:
      break;
      
    case Primitive:
      result = primitiveToConservative(result);
      break;
      
    case Characteristic:
      result = rightEigenMatrix(averageC,axisNo)*result;
  };
  
  return result;
}

/*
_______________________________________________________________________________

Compute adiabatic wall condition
_______________________________________________________________________________

*/

void Cell::setAdiabaticWallAverage(const Vector& average) 
{
  int i;
  Vector V(Dimension);
  real p, T, rho, rhoE, rhoY;

  rho = average.at(1);
  
  for (i=1; i<= Dimension; i++)
    V.set(i,average.at(i+1)/rho);
    
  rhoE = average.at(Dimension+2);
  rhoY = (SpeciesNb > 1)? average.at(Dimension+3):0.;
  
  p = EOS_P(rho, V, rhoE, rhoY);
  T = MolarMass*p/(rho*R);
  
  // --- Compute new values ---
  
  V = -V;
  rho = MolarMass*p/(R*T);
  rhoE = EOS_rhoE(rho, V, p, rhoY);
  
  // --- set new values ---
  
  setAverage(1,rho);
  
  for (i=1; i<= Dimension; i++)
    setAverage(i+1,rho*V.at(i));
  
  setAverage(Dimension+2, rhoE);

  if (SpeciesNb > 1)
  {
    for (i=1; i<=SpeciesNb-1; i++)
      setAverage(Dimension+2+i, average.at(Dimension+2+i));
  }

  
}  

/*
_______________________________________________________________________________

Compute isothermal wall condition
_______________________________________________________________________________

*/

void Cell::setIsothermalWallAverage(const Vector& average) 
{
  int i;
  Vector V(Dimension);
  real p,Tw, rho, rhoE, rhoY;

  rho = average.at(1);
  
  for (i=1; i<= Dimension; i++)
    V.set(i,average.at(i+1)/rho);
    
  rhoE = average.at(Dimension+2);
  rhoY = (SpeciesNb > 1)? average.at(Dimension+3):0.;
  
  p = EOS_P(rho, V, rhoE, rhoY);
  
  // --- Choose the constant wall temperature ---
  
  Tw = T();
  
  // --- Compute new values ---
  
  rho = MolarMass*p/(R*Tw);
  rhoE = EOS_rhoE(rho, V, p, rhoY);
  
  // --- set new values ---
  
  setAverage(1,rho);
  
  for (i=1; i<= Dimension; i++)
    setAverage(i+1,-rho*V.at(i));
  
  setAverage(Dimension+2, rhoE);

  if (SpeciesNb > 1)
  {
    for (i=1; i<=SpeciesNb-1; i++)
      setAverage(Dimension+2+i, average.at(Dimension+2+i));
  }

  
}  

/*
_______________________________________________________________________________

Compute inlet condition
_______________________________________________________________________________

*/

void Cell::setInletAverage(const Vector& average) 
{
  int i;
  Vector V(Dimension);
  real p, rho, rhoE, rhoY, T;

  // --- Extract pressure from average
  
  rho = average.at(1);
  
  for (i=1; i<= Dimension; i++)
    V.set(i,average.at(i+1)/rho);
      
  rhoE = average.at(Dimension+2);
  rhoY = (SpeciesNb > 1)? average.at(Dimension+3):0.;

  p = EOS_P(rho, V, rhoE, rhoY);
  
  // --- Get temperature and velocity from the boundary
 
  T = temperature();
  V = velocity();
  
  // --- In case of supersonic flow, get also the pressure
  
  if (Regime == Supersonic) 
    p = pressure();
  
  // --- Compute density ---
  
  rho = MolarMass*p/(R*T);
  
  // --- Store new energy in this cell ---
  
  rhoE = EOS_rhoE(rho, V, p, rhoY);
  
  setAverage(1, 1, rho);
  
  for (i=1; i<= Dimension; i++)
    setAverage(i+1,1,rho*V.at(i));
  
  setAverage(Dimension+2, 1, rhoE);
  if (SpeciesNb > 1) setAverage(Dimension+3,1,rhoY);
}  

/*
_______________________________________________________________________________

Compute outlet condition
_______________________________________________________________________________

*/

void Cell::setOutletAverage(const Vector& average1, const Vector& average2) 
{
  int i;
  Vector V1(Dimension), V2(Dimension), V(Dimension);
  real rho1, rhoE1, rhoY1;
  real rho2, rhoE2, rhoY2;
  real p, rho, rhoE, rhoY, T;


  // --- Extract first state from cell1 ---
  
  rho1 = average1.at(1);
  
  for (i=1; i<= Dimension; i++)
    V1.set(i,average1.at(i+1)/rho1);
      
  rhoE1 = average1.at(Dimension+2);
  rhoY1 = (SpeciesNb > 1)? average1.at(Dimension+3):0.;

  // --- Extract second state from cell2 ---

  rho2 = average2.at(1);
  
  for (i=1; i<= Dimension; i++)
    V2.set(i,average2.at(i+1)/rho2);
      
  rhoE2 = average2.at(Dimension+2);
  rhoY2 = (SpeciesNb > 1)? average2.at(Dimension+3):0.;
  
  // Choose the second state
  
  p = EOS_P(rho2, V2, rhoE2, rhoY2);
  rho = rho2;
  rhoY = rhoY2;
  V = V2;
  
  // --- Store new energy in this cell ---
 
  rhoE = EOS_rhoE(rho, V, p, rhoY);
  setAverage(1, 1, rho);
  
  for (i=1; i<= Dimension; i++)
    setAverage(i+1,1,rho*V.at(i));
  
  setAverage(Dimension+2, 1, rhoE);
  if (SpeciesNb > 1) setAverage(Dimension+3,1,rhoY);
}  

/*
_______________________________________________________________________________

Compute outlet condition
_______________________________________________________________________________

*/

void Cell::setFreeSlipAverage(const Vector& average, const int axisNo) 
{
  int i;
  Vector V(Dimension);
  real p, rho, rhoE, rhoY, T;

  // --- Extract pressure from average ---
  
  rho = average.at(1);
  
  for (i=1; i<= Dimension; i++)
    V.set(i,average.at(i+1)/rho);
      
  rhoE = average.at(Dimension+2);
  rhoY = (SpeciesNb > 1)? average.at(Dimension+3):0.;

  p = EOS_P(rho, V, rhoE, rhoY);
  
  T = MolarMass*p/(rho*R);
  
  // --- Invert the component of the velocity in the direction axisNo ---
  
  V.set(axisNo,-V.at(axisNo));
  
  // --- Compute density ---
  
  rho = MolarMass*p/(R*T);  
    
  // --- Store new energy in this cell ---
 
  rhoE = EOS_rhoE(rho, V, p, rhoY);
  setAverage(1, 1, rho);
  
  for (i=1; i<= Dimension; i++)
    setAverage(i+1,1,rho*V.at(i));
  
  setAverage(Dimension+2, 1, rhoE);
  if (SpeciesNb > 1) setAverage(Dimension+3,1,rhoY);
}  




/*
_______________________________________________________________________________

Compute source term with 6th order Runge-Kutta method (for detonation)
_______________________________________________________________________________

*/

void Cell::computeSource()
{ 
  // --- Init ---
  
  const real c1 = 5179.0/57600.0;
  const real c3 = 7571.0/16695.0;
  const real c4 = 393.0/640.0;
  const real c5 = -92097.0/339200.0;
  const real c6 = 187.0/2100.0;
  const real c7 = 1.0/40.0;
  
  const real d1 = 35.0/384.0;
  const real d3 = 500.0/1113.0;
  const real d4 = 125.0/192.0;
  const real d5 = -2187.0/6784.0;
  const real d6 = 11.0/84.0;

  const real b21 = 1.0/5.0;
  const real b31 = 3.0/40.0;
  const real b32 = 9.0/40.0;
  const real b41 = 44.0/45.0;
  const real b42 = -56.0/15.0;
  const real b43 = 32.0/9.0;
  const real b51 = 19372.0/6561.0;
  const real b52 = -25360.0/2187.0;
  const real b53 = 64448.0/6561.0;
  const real b54 = -212.0/729.0;
  const real b61 = 9017.0/3168.0;
  const real b62 = -355.0/33.0;
  const real b63 = 46732.0/5247.0;
  const real b64 = 49.0/176.0;
  const real b65 = -5103.0/18656.0;
  const real b71 = 35.0/384.0;
  const real b72 = 0.0;
  const real b73 = 500.0/1113.0;
  const real b74 = 125.0/192.0;
  const real b75 = -2187.0/6784.0;
  const real b76 = 11.0/84.0;
  
  real dt, rt; 			// local time step, remaining time
  real rhoY, rhoY_low; 		// 5th order solution of rhoY, and 4th order
  real error = 0.;
  real tolerance = RKFTolerance;
  real s;
  
  real k1, k2, k3, k4, k5, k6, k7;
  real y, y0;

  real a = MolarMass*(Gamma-1.)/R;
  real b, c, Te;
  
  int i=0;

  // --- Execution ---

  rhoY = average(Dimension+3,1);
  b = Q0/rho();
  c = e()-0.5*velocity()*velocity();
  dt = RKFTimeStepFactor*TimeStep;
  rt = 0.5*TimeStep;      
  
  while (rt > 0.)
  {
    i++;
    rt -= dt;
    y0   = rhoY;
       
    // Explicit 4-5th order Runge-Kutta method
    
    y = y0;
    Te = a*(c-b*y);
    k1 = (Te > Ti) ? -K0*y*exp(-Ta/Te):0.;
    
//    if (i == 1) printf("x = %15.8e, rhoY = %15.8e, T = %15.8e, k1 = %15.8e\n",x(), y, Te, k1);
    
    y = y0 + dt * b21* k1;
    Te = a*(c-b*y);
    k2 = (Te > Ti) ? -K0*y*exp(-Ta/Te):0.;
    
    y = y0 + dt * ( b31*k1 + b32*k2);
    Te = a*(c-b*y);
    k3 = (Te > Ti) ? -K0*y*exp(-Ta/Te):0.;
   
    y = y0 + dt * ( b41*k1 + b42*k2 + b43*k3);
    Te = a*(c-b*y);
    k4 = (Te > Ti) ? -K0*y*exp(-Ta/Te):0.;
    
    y = y0 + dt * ( b51*k1 + b52*k2 + b53*k3 + b54*k4);
    Te = a*(c-b*y);
    k5 = (Te > Ti) ? -K0*y*exp(-Ta/Te):0.;
    
    y = y0 + dt * ( b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5);
    Te = a*(c-b*y);
    k6 = (Te > Ti) ? -K0*y*exp(-Ta/Te):0.;
    
    y = y0 + dt * ( b71*k1 + b72*k2 + b73*k3 + b74*k4 + b75*k5 + b76*k6);
    Te = a*(c-b*y);
    k7 = (Te > Ti) ? -K0*y*exp(-Ta/Te):0.;
      
    rhoY_low  = y0 +  dt * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6 + c7 * k7);
    rhoY      = y0 +  dt * (d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6);

    error = fabs(rhoY-rhoY_low)/fabs(rhoY);
    if (rhoY == 0. && rhoY_low == 0.) error = 0.;

    // Runge-Kutta-Fehlberg time step computation

    s = exp(0.25*log((tolerance*dt)/(error*TimeStep)));
    if (error == 0.) s = 3.;
    
    // Adapt time step
    
    if (s > 2) dt *= 2.;
    if (s < 1) dt *= 0.5;
    if (dt > rt) dt = rt;
  }
  
  if (rhoY > rho()) rhoY = rho();
  if (rhoY < 0.) rhoY = 0.;
  
  setAverage(Dimension+3, 1, rhoY);
  
}

/*
_______________________________________________________________________________

Penalize momentum
_______________________________________________________________________________

*/

void Cell::penalizeMomentum()
{
  int i;
  Vector M(Dimension);
  real k = exp(-0.25*TimeStep/Permeability);
  
#ifdef debug_mode  
  if (Equation != NavierStokes)
  {
    cout << "lara: in function Cell::penalize()" << endl;
    cout << "lara: equation must be Navier-Stokes" << endl;
    exit(-1);
  }
#endif
    
    
   M.set(1,density()*PenalVelocity);
   for (i=1;i<=Dimension;i++)
     setMomentum(i,M.at(i) + (momentum(i) - M.at(i))*k);
}

/*
_______________________________________________________________________________

Penalize energy
_______________________________________________________________________________

*/

void Cell::penalizeEnergy()
{
  real rhoE;
  real k = exp(-0.25*TimeStep/Permeability);
      
   rhoE = density()*R*PenalTemperature/((Gamma-1.)*MolarMass) + 0.5*density()*SQR(PenalVelocity);
   setAverage(Dimension+2, 1, rhoE + (average(Dimension+2,1) - rhoE)*k);  
   
}

/*
_______________________________________________________________________________

Penalize density
_______________________________________________________________________________

*/

void Cell::penalizeDensity(const Cell* CIm3, const Cell* CIm2, const Cell* CIm1, const Cell* CIp1, const Cell* CIp2, const Cell* CIp3, const int axisNo)
{
//	real fluxIp12; // f(i+1/2)
//	real fluxIm12; // f(i-1/2)
	real Dt = 0.5*TimeStep;
	real Dx = size(axisNo);
	
//	fluxIp12 = fluxRoeDensity(CIm2,CIm1,this,CIp1,CIp2,CIp3,axisNo);
//	fluxIm12 = fluxRoeDensity(CIm3,CIm2,CIm1,this,CIp1,CIp2,axisNo);
	
//	setDensity(density()-(1./Porosity-1.)*Dt/Dx*(fluxIp12-fluxIm12));
	
	setDensity(density()-(1./Porosity-1.)*Dt/(2.*Dx)*(CIp1->momentum(axisNo)-CIm1->momentum(axisNo)));
	
}

/*
_______________________________________________________________________________

Conservative to primitive
_______________________________________________________________________________

*/

Vector Cell::conservativeToPrimitive(const Vector& average) const
{
  int i;
  real rho;
  Vector V;
  real rhoE;
  real p;
  Vector result(QuantityNb);
  
  switch(Equation)
  {
    case ConvectionDiffusion:
    case Burgers:
    default:
      result = average;
      break;
      
    case NavierStokes:
      
      // Get conservative variables and compute primitive variables
	
      rho = average.at(1);
	
      for (i=1; i<=Dimension;i++)
	V.set(i,average.at(i+1)/rho);
	
      rhoE = average.at(Dimension+2);
      p = (Gamma-1)*(rhoE - 0.5*rho*V*V);
	
      // Put primitive variables
	
      result.set(1,rho);
	
      for (i=1; i<=Dimension; i++)
	result.set(i+1,V.at(i));
	
      result.set(Dimension+2,p);
      
      break;
  };
  
  return result;
  
}

/*
_______________________________________________________________________________

Primitive to conservative
_______________________________________________________________________________

*/

Vector Cell::primitiveToConservative(const Vector& average) const
{
  int i;
  real rho;
  Vector V;
  real rhoE;
  real p;
  Vector result(QuantityNb);
  
  switch(Equation)
  {
    case ConvectionDiffusion:
    case Burgers:
    default:
      result = average;
      break;
      
    case NavierStokes:
      
      // Get conservative variables and compute primitive variables
	
      rho = average.at(1);
	
      for (i=1; i<=Dimension;i++)
	V.set(i,average.at(i+1));
	
      p = average.at(Dimension+2);
      rhoE = p/(Gamma-1) + 0.5*rho*V*V;
	
      // Put primitive variables
	
      result.set(1,rho);
	
      for (i=1; i<=Dimension; i++)
	result.set(i+1,rho*V.at(i));
	
      result.set(Dimension+2,rhoE);
      
      break;
  };
  
  return result;
  
}

/*
_______________________________________________________________________________

Link boundary cell (l,i1,j1,k1) with cell (l, i2, j2, k2)
_______________________________________________________________________________

*/

void Cell::linkBoundaryCell(const int BC, Cell* cell1, Cell* cell2, Cell* cell3, const int axisNo)
{
  switch (BC)
  {
    case Inlet:
    default:
      setInletAverage(cell1->average());
      break;
    
    case Outlet:
      setOutletAverage(cell1->average(),cell2->average(StorageNb));
      break;
      
    case Periodic:
      setAverage(cell3->average());
      break;
    
    case FreeSlip:
      if (Equation==NavierStokes) 
	setFreeSlipAverage(cell1->average(),axisNo);
      else
	setAverage(cell1->average());
      break;
      
    case AdiabaticWall:
      setAdiabaticWallAverage(cell1->average());
      break;
      
    case IsothermalWall:
      setIsothermalWallAverage(cell1->average());
      break;
  };
}
