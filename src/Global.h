// --- Precision ---

#ifdef double_precision
	#define real double
#else
	#define real float
#endif

// --- FUNCTIONS ---

#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))
#define SGN(x) (((x)>0.)? 1.:(((x)<0.)?-1.:0.))
#define ABS(x) (((x)>0.)? (x):-(x))

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)<(y))?(y):(x))

// EQUATION
#define NavierStokes 1
#define ConvectionDiffusion 2
#define Burgers 3

// BOUNDARY CONDITIONS
#define Periodic 1
#define FreeSlip 2
#define Inlet 3
#define Outlet 4
#define AdiabaticWall 5
#define IsothermalWall 6

// SCHEMES
#define McCormack 1
#define Roe 2
#define AUSM 3

// RECONSTRUCTION
// None = 1
#define MUSCL 2
#define WENO 3

// LIMITERS
#define None 1
#define MinMod 2
#define VanAlbada 3
#define VanLeer 4
#define Ospre 5
#define Koren 6
#define SuperBee 7

// DECENTERING OPTIONS
#define Decentered 1
#define Centered 2
#define ThirdOrder 3

// FILE FORMAT OPTIONS
#define Gnuplot 1
#define VTK 2

// FILE TYPE OPTIONS
#define ASCII 1
#define Binary 2

// TREE OPTIONS
// None = 1
#define Node 2
#define Leaf 3
#define Boundary 4

// REGIME OPTIONS
#define Supersonic 1
#define Subsonic 2

// WENO VARIABLES
#define Conservative 1
#define Primitive 2
#define Characteristic 3

