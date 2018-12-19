class Matrix
{
public:
	// Constructors and distructor

	Matrix();
	Matrix(const int n, const int m);
	Matrix(const Matrix& M);
	~Matrix();

	// Set functions

	inline void set(const int i, const int j, const real x)
	{ _U[(i-1)*_columns+(j-1)] = x; }

	void operator= (const Matrix& M);

	void resize(const int n, const int m);

	void raise(const int n, const int m);
	
	void reset();

	// Get functions

	inline real at(const int i, const int j) const
	{ return _U[(i-1)*_columns+(j-1)]; } 

	inline int columns() const
	{ return _columns; }

	inline int lines() const
	{ return _lines; }

	Vector vectorColumn(const int column) const;

	// Additions

	void operator+= (const Matrix& M);
	Matrix operator+ (const Matrix& M) const;

	// Substractions

	void operator-= (const Matrix& M);
	Matrix operator- (const Matrix& M) const;

	// Multiplication with real

	void operator*= (const real x);

	// Internal product

	Matrix operator* (const Matrix& M) const;

	// Matrix-Vector product

	Vector operator* (const Vector& V) const;

	// swap lines
	
	void swapLines(const int line1, const int line2);
	
	// swap columns
	
	void swapColumns(const int column1, const int column2);
	
private:
	real* _U;
	int _lines, _columns;
};

/*
_______________________________________________________________________________

External functions
_______________________________________________________________________________

*/

Matrix operator* (const real x, const Matrix& M);

Matrix fabs(const Matrix& M);