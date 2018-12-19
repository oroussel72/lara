/*
_______________________________________________________________________________

List of parameters
_______________________________________________________________________________

Geometry
_______________________________________________________________________________

*/

extern int Dimension;	// Dimension of the problem
extern int ScaleNb;	// Number of scales (resolution = 2^(ScaleNb*Dimension))

extern int NX;
extern int NY;	// Number of cells in each direction for the FV solver
extern int NZ;


extern real XMin, XMax;	// Boundaries in x
extern real YMin, YMax;	// Boundaries in y
extern real ZMin, ZMax;	// Boundaries in z

extern int BCXMin, BCXMax; // Boundary conditions in XMin, XMax
extern int BCYMin, BCYMax; // Boundary conditions in YMin, YMax
extern int BCZMin, BCZMax; // Boundary conditions in ZMin, ZMax

// Accepted values: Periodic, FreeSlip, NoSlip, Inlet, Outlet

/*
_______________________________________________________________________________

Multiresolution
_______________________________________________________________________________

*/

extern bool Multiresolution; // true = MR, false = FG
extern real Tolerance;
extern int ScaleMin;	// Minimal scale where computations are performed
extern bool RefineInCenter; 	// true = keep refined in the center (for Sedov)
extern int RefineWidth;		// width of center refinement
extern bool UseVariableTolerance;	// yes = use variable tolerance

/*
_______________________________________________________________________________

Numerics
_______________________________________________________________________________

*/

extern int Equation;   // ConvectionDiffusion, Burgers or NavierStokes
extern int Regime; 	// Supersonic or Subsonic (for the inlet and outlet BC, Navier-Stokes only)
extern int Scheme;	// McCormack, Roe, AUSM
extern int Reconstruction; // MUSCL, WENO
extern int Limiter;	// VanAlbada, VanLeer, Ospre, Koren (MUSCL only)
extern int Centering; 	// Centered, Decentered, ThirdOrder (MUSCL only)
extern int WENOVariables; // Conservative, Primitive, Characteristic (WENO only)
extern int QuantityNb;	// Number of conservative quantities
extern int SpeciesNb;	// Number of chemical species	
extern int StorageNb;	// Number of storages for the cell-average values
extern int StepNb;	// Number of steps for a multi-stage time integration
extern real TimeStep;	// Time step
extern real CFL;	// Courant-Friedrich-Lewy number
extern real TotalTime;	// Total time of the computation
extern int WriteNb;	// Number of data to write in file
extern int WriteFormat;	// Gnuplot, VTK
extern int WriteType;	// ASCII or Binary
extern bool WriteStrip; // Write date into strip format
extern bool WriteStripCoarse; // Unrefine the strip format from 1 level (only for FG computations)
extern bool Restart;	// true = restart from backup date, false = start from initial condition

/*
_______________________________________________________________________________

Penalization
_______________________________________________________________________________

*/

extern bool UsePenalization;
extern real Permeability;		// eta coefficient of penalization
extern real Porosity;
extern real PenalTemperature;
extern real PenalVelocity;

/*
_______________________________________________________________________________

Physics
_______________________________________________________________________________

*/

extern bool IsInviscid;		// true = inviscid flow
extern real Diffusivity;	// diffusivity for ConvectionDiffusion/Burgers
extern real Celerity;		// celerity for Convection/Diffusion
extern bool ConstantViscosity; // false = use Sutherland's law
extern real Viscosity;
extern Vector MassDiffusivity;
extern real Gamma;		// isentropic exponent
extern real R;			// ideal gas constant
extern real MolarMass;
extern real Pr;			// Prandtl number
extern real Tr;			// Reference temperature (Sutherland's law)
extern real Ts;			// Sutherland temperature (Sutherland's law)

/*
_______________________________________________________________________________

Detonation
_______________________________________________________________________________

*/

extern bool Detonation;		// true = use detonation
extern real Ta;			// Activation temperature
extern real Ti;			// Ignition temperature
extern real Q0;			// Heat release
extern real K0;			// Rate constant
extern real RKFTolerance;	// Tolerance for the RKF procedure
extern real RKFTimeStepFactor; // Time step factor for the RKF procedure

/*
_______________________________________________________________________________

Internal parameters
_______________________________________________________________________________

*/

extern int RefreshNb;		// Number of screen refreshments
extern int IterationNo;		// current iteration number
extern int StepNo;		// step number in multi-stage time integrations
extern real AverageTimeStep;   // average time step
extern real RemainingTime;	// remaining time
extern int WriteNo;		// current write data number
extern bool Verbose;		// true = verbose mode
extern real SpaceStep;		// smallest space step

extern Timer LaraTimer;		// Timer for the code
extern real Conductivity;	// Thermal conductivity
extern Vector MaxAverage;	// Maximal values of the cell-averages
extern int NodeNb;		// Number of nodes and leaves
extern int LeafNb;		// Number of leaves
extern real AverageNodeNb;	// Average number of nodes and leaves
extern real AverageLeafNb;	// Average number of leaves

extern Vector PredictionError;		// prediction error
extern Vector PredictionErrorOld;	// old prediction error
extern real PredictionFactor;		// power to affect to prediction error change
extern real InitialTolerance;		// Initial tolerance
extern real MultiresolutionError;	// Multiresolution error
extern real StartTime;			// Time where the computation starts

/*
_______________________________________________________________________________

Integral values
_______________________________________________________________________________

*/

extern real MaxEigenvalue;		// Maximal eigenvalue
extern real KineticEnergy;		// Kinetic energy
extern real MaxDiffusivity; 		// Maximal diffusivity 
extern real TimeAverageMaxDensity;	// Maximal density (averaged in time)


/*
_______________________________________________________________________________

Routines
_______________________________________________________________________________

*/

void initParameters();

// Initial condition function

Vector initCondition(const real x, const real y, const real z);

// True is the point (x,y,z) is immersed in the boundary

bool isImmersed(const real x, const real y, const real z, const real t);

// Return true when it is time to show evolution on screen

bool timeToShow();

// Return true when it is time to write output data

bool timeToWrite();

// Show the performance of the computation on screen or on file performance.out

void showPerformance();
