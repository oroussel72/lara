class RegularMesh
{
public:
	// Constructor and distructor

	RegularMesh();
	~RegularMesh();

	// Pointer to cell
	
	inline Cell* cell(const int i, const int j, const int k) const
	{ return (_cell + i + (NX+6)*j + SQR(NX*NY+6)*k);}

	// Write data
	
	inline void writeLn(FILE* f, const real arg) const
	{ 
	  double x = arg;
	  if (WriteType == ASCII)
	    fprintf(f,"%23.16e\n",x);
	  else
	    fwrite(&x,sizeof(double),1,f);
	}
	
	inline void write(FILE* f, const real arg) const
	{ 
	  double x = arg;
	  if (WriteType == ASCII)
	    fprintf(f,"%23.16e  ",x);
	  else
	    fwrite(&x,sizeof(double),1,f);
	}
	// Compute the divergence term, i.e. the fluxes in every direction
	
	void computeDivergence();
	
	// Compute the gradients in every direction
	
	void computeGradient();
	
	// Compute the source term
	
	void computeSource();
	
	// Copy the cell-average values from source to target
	
	void copyAverage(const int sourceStorageNo, const int targetStorageNo);

	// Store cell-average values from stage 1 to stage 2 for multi-stage computations
	
	inline void storeAverage()
	{ copyAverage(1,2); }
	
	inline void storeOldAverage()
	{ copyAverage(StorageNb-1,StorageNb); }
	
	// Backup cell-average values
	
	void backup();
	
	// Restore cell-average values
	
	void restore();
	
	// Penalize cells in solid
	
	void penalize();
	
	// Interpolate fluid when going out of the solid
	
	void interpolateFluid(const int i, const int j, const int k);
	
	// Interpolate solid when going out of the fluid
	
	void interpolateSolid(const int i, const int j, const int k);
	
	// Perfom all these interpolations (fluid and solid)
	
	void render();
	
	// Compute one step of the multi-stage time integration method
	
	void computeStep();

	// Compute integral values
	
	void computeIntegralValues();
	
	// Compute time step
	
	void computeTimeStep();
	
	// Check the stability of the computation
	
	void checkStability();
	
	// Set the boundary conditions
	
	void setBoundaryConditions();
	
	// Compute complete iteration
	
	void iterate();
	
	// Write cell-average values into data file

	void writeAverage();
			
	// Write integral values into data file

	void writeIntegral();
	
	// Compute initial condition
	
	void initAverage();
	
	// Init mesh
	
	void init();
	
	// End computation
	
	void end();
	
	
private:
	Cell* _cell;			// Array of cells
};
