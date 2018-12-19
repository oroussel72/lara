class Cell
{
public:
	// Constructor and distructor

	Cell();
	~Cell();
	
	// Set procedures
	
	inline void setCenter(const int axisNo, const real x)
	{ _center.set(axisNo,x); }

	inline void setSize(const int axisNo, const real x)
	{ _size.set(axisNo,x); }

	inline void setAverage(const int quantityNo, const int storageNo, const real x)
	{ _average.set(quantityNo,storageNo,x); }

	inline void setAverage(const int quantityNo, const real x)
	{ _average.set(quantityNo,1,x); }	

	inline void setDensity(const real x)
	{ _average.set(1,1,x); }
	
	inline void setMomentum(const int axisNo, const real x)
	{ _average.set(axisNo+1,1,x); }
	
	inline real EOS_P(const real rho, const Vector& V, const real rhoE, const real rhoY) const
	{ return (Gamma-1.)*(rhoE - 0.5*rho*V*V - ((Detonation)? rhoY*Q0:0.));}
	
	inline real EOS_rhoE(const real rho, const Vector& V, const real p, const real rhoY) const
	{ return p/(Gamma-1.) + 0.5*rho*V*V + ((Detonation)? rhoY*Q0:0.);}
	
	void setAverage(const int storageNo, const Vector& V);

	inline void setAverage(const Vector& V)
	{ setAverage(1,V); }
	
	void resetAverage(const int storageNo=1);
	
	inline void setGradient(const int quantityNo, const int axisNo, const real x)
	{ _gradient.set(quantityNo,axisNo,x); }

	inline void setDivergence(const Vector& V)
	{ _divergence = V;}

	inline void addDivergence(const Vector& V)
	{ _divergence += V;}
	
	inline void resetDivergence()
	{ _divergence.reset();}

	// Penalization procedures
	
	inline bool wasInSolid() const
	{ return isImmersed(x(),y(),z(), TotalTime-RemainingTime-TimeStep); };
	
	inline bool isInSolid() const
	{ return isImmersed(x(),y(),z(), TotalTime-RemainingTime);}
	
	inline bool isInSolidAtMid() const
	{ return isImmersed(x(),y(),z(), TotalTime-RemainingTime + 0.5*TimeStep);}

	inline bool isInSolidAtQuarter() const
	{ return isImmersed(x(),y(),z(), TotalTime-RemainingTime + 0.25*TimeStep);}

	inline bool isInSolidAtThreeQuarter() const
	{ return isImmersed(x(),y(),z(), TotalTime-RemainingTime + 0.75*TimeStep);}

	void penalizeMomentum();

	void penalizeDensity(const Cell* CIm3, const Cell* CIm2, const Cell* CIm1, const Cell* CIp1,const Cell* CIp2, const Cell* CIp3,const int axisNo);

	void penalizeEnergy();

	// Get procedures

	inline real center(const int axisNo) const
	{ return _center.at(axisNo); }

	inline Vector center() const
	{ return _center; }

	inline real x() const
	{ return _center.at(1); }

	inline real y() const
	{ return (Dimension > 1) ? _center.at(2):0.; }

	inline real z() const
	{ return (Dimension > 2) ? _center.at(3):0.; }

	inline real size(const int axisNo) const
	{ return _size.at(axisNo); }

	inline Vector size() const
	{ return _size; }

	inline real dx() const
	{ return _size.at(1); }

	inline real dy() const
	{ return (Dimension > 1) ? _size.at(2):1.; }

	inline real dz() const
	{ return (Dimension > 2) ? _size.at(3):1.; }

	inline real average(const int quantityNo, const int storageNo) const
	{ return _average.at(quantityNo, storageNo); }

	inline Vector average(const int storageNo=1) const
	{ return _average.vectorColumn(storageNo); }
	
	inline real gradient(const int quantityNo, const int axisNo) const
	{ return _gradient.at(quantityNo, axisNo);}

	inline Matrix gradient() const
	{ return _gradient;}
	
	inline Vector divergence() const
	{ return _divergence;}

	inline real density(const int storageNo=1) const
	{ return _average.at(1, storageNo);}

	inline real rho() const
	{ return density(1); }
	
	inline real energy(const int storageNo=1) const
	{ return _average.at(Dimension+2, storageNo);} 

	inline real e() const
	{ return _average.at(Dimension+2,1)/_average.at(1,1);}
	
	inline real velocity(const int axisNo, const int storageNo) const
	{ return (Equation != NavierStokes)? _average.at(axisNo,storageNo):_average.at(axisNo+1, storageNo)/_average.at(1,storageNo);}

	inline real momentum(const int axisNo, const int storageNo=1) const
	{ return _average.at(axisNo+1,storageNo);}
	
	real momentumNorm(const int storageNo=1) const;
	
	inline real u() const
	{ return velocity(1,1);}
	
	inline real v(const int axisNo=2) const
	{ return velocity(axisNo,1); }
	
	inline real w() const
	{ return velocity(3,1);}
	
	Vector velocity(const int storageNo=1) const;
	
	real velocityNorm(const int storageNo=1) const; 

	inline real partialMass(const int speciesNo, const int storageNo=1) const
	{ return _average.at(2+Dimension+speciesNo, storageNo)/density(storageNo);}

	inline real Y(const int speciesNo) const
	{ return partialMass(speciesNo,1); }
	
	real temperature(const int storageNo=1) const;

	inline real T() const
	{ return temperature(1); }
	
	bool isNaN() const;
	
	real pressure(const int storageNo=1) const;
	
	inline real p() const
	{ return pressure(1); }
	
	inline real speedOfSound(const int storageNo=1) const
	{ return sqrt(Gamma*pressure(storageNo)/density(storageNo));}

	inline real c() const
	{ return speedOfSound(1);}

	inline real sutherland(const real T) const
	{ return (ConstantViscosity)? Viscosity:Viscosity*exp(1.5*log(T/Tr))*(Tr+Ts)/(T+Ts);}
	
	real maxEigenvalue(const int storageNo=1) const;
	
	Vector eulerFlux(const Vector& average, const int axisNo) const; 

	void computeSource();
	
	Matrix leftEigenMatrix(const Vector& average, const int axisNo) const;
	
	Matrix rightEigenMatrix(const Vector& average, const int axisNo) const;
	
	Matrix absDiagonalMatrix(const Vector& average, const int axisNo) const;
	
	// Copy procedures

	void copyStorage(const int sourceStorageNo, const int targetStorageNo);

	// Compute gradient of stage 1

	void computeGradient(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const int axisNo);

	// Compute flux

	Vector inviscidFlux(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo) const;

	Vector viscousFlux(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo) const;

	// Numerical schemes

	Vector schemeMcCormack(const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const int axisNo) const;

	Vector schemeRoe(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo) const;
	
	real schemeRoeDensity(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo) const;

	Vector schemeAUSM(const Cell* Cim2, const Cell* Cim1, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo) const;
	
	// Reconstruction procedures
	
	Vector roeAverage(const Vector& average1, const Vector& average2) const;

	Vector leftLimitedAverage(const Vector& averageIm1, const Vector& averageI, const Vector& averageIp1) const;
	
	Vector rightLimitedAverage(const Vector& averageI, const Vector& averageIp1, const Vector& averageIp2) const;
	
	real limiter(const real r) const;
	
	Vector WENOAverage(const Vector& averageIm2, const Vector& averageIm1, const Vector& averageI, const Vector& averageIp1, const Vector& averageIp2, const int axisNo) const;

	Vector conservativeToPrimitive(const Vector& average) const;
	
	Vector primitiveToConservative(const Vector& average) const;
			
	// Tree procedures
	
	inline void setStatus(const int value)
	{ _status = value;}
	
	inline int status() const
	{ return _status;}
	
	// Boundary procedures
	
	void linkBoundaryCell(const int BC, Cell* cell1, Cell* cell2, Cell* cell3, const int axisNo);
	
	void setAdiabaticWallAverage(const Vector& average);
	
	void setIsothermalWallAverage(const Vector& average);

	void setInletAverage(const Vector& average);
	
	void setOutletAverage(const Vector& average1, const Vector& average2);
	
	void setFreeSlipAverage(const Vector& average, const int axisNo);
	
	
private:
	Vector _center;		// coordinates of the center of the cell
	Vector _size;		// size of the cell (in x, y, z)
	Vector _divergence;	// divergence vector
	Matrix _average;	// vector of the cell-average values in different storages
	Matrix _gradient;	// velocity gradient
	int _status;		// None, Node or Leaf
	
};
 
/*
_______________________________________________________________________________

EXTERNAL FUNCTIONS
_______________________________________________________________________________

*/

inline Vector flux(const Cell* Cim2, const Cell* Cim1, const Cell* Ci, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo)
{ return (IsInviscid) ? Ci->inviscidFlux(Cim2, Cim1, Cip1, Cip2, Cip3, axisNo) : 
  Ci->inviscidFlux(Cim2, Cim1, Cip1, Cip2, Cip3, axisNo) + Ci->viscousFlux(Cim2, Cim1, Cip1, Cip2, Cip3, axisNo);} 

inline real fluxRoeDensity(const Cell* Cim2, const Cell* Cim1, const Cell* Ci, const Cell* Cip1, const Cell* Cip2, const Cell* Cip3, const int axisNo)
{ return Ci->schemeRoeDensity(Cim2, Cim1, Cip1, Cip2, Cip3, axisNo);}