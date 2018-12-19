class AdaptiveMesh
{

public:

  // Constructor and distructor
  
  AdaptiveMesh();
  ~AdaptiveMesh();	

  // Number of cells on each level
  
  inline int N(const int l) const
  { return (1<<l); }
    
  // Returns the cells at level l and position (i,j,k)
  
  inline Cell* cell(const int l, const int i, const int j, const int k) const
  { return (_cell + _cellNb[l] + (i+3) + ((Dimension > 1) ? (N(l)+6)*(j+3):0) + ((Dimension > 2) ? SQR((N(l))+6)*(k+3):0));}
  
  inline Cell* parent(const int l, const int i, const int j, const int k) const
  { return cell(l-1,i/2,j/2,k/2);}
  
  inline Cell* uncle(const int l, const int i, const int j, const int k, const int ni, const int nj, const int nk) const
  { return cell(l-1,i/2+ni,j/2+nj,k/2+nk);}  
  
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
  // Local procedures
    
  void split(const int l, const int i, const int j, const int k);
  
  void combine(const int l, const int i, const int j, const int k);
  
  Vector predictCell(const int l, const int i, const int j, const int k) const;
  
  void projectCell(const int l, const int i, const int j, const int k);
 
  Vector detail(const int l, const int i, const int j, const int k) const;

  bool isUncleOfNode(const int l, const int i, const int j, const int k) const;
  
  bool isCousinOfLeaf(const int l, const int i, const int j, const int k) const;

  bool isParentOfLeaf(const int l, const int i, const int j, const int k) const;

  bool isSmall(const int l, const Vector& average) const;
  
  void interpolateFluid(const int l, const int i, const int j, const int k);

  void interpolateSolid(const int l, const int i, const int j, const int k);

  // Global procedures
  
  void copyAverage(const int sourceStorageNo, const int targetStorageNo);
  
  // Store cell-average values from stage 1 to stage 2 for multi-stage computations
  
  inline void storeAverage()
  { copyAverage(1,2); }

  void countLeavesAndCells();
  
  void setBoundaryConditions();
  
  void checkStability();
  
  void update();
  
  void updateAll();

  void adapt();
  
  void computeDivergence();
  
  void computeStep();
  
  void computeIntegralValues();
  
  void computeTimeStep();
  
  void computeGradient();
  
  void computeSource();
  
  void writeAverage();
  
  void writeTree();
  
  void writeIntegral();
  
  void writeStrip();
  
  void backup();
  
  void restore();
  
  void penalize();
  
  void render();
  
  void initAverage();
  
  void init();
  
  void end();
  
  void iterate();

private:
  
  Cell* _cell;	// array of cells
  int* _cellNb; // number of cells from level 0 to level l-1
};
