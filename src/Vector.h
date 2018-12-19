class Vector
{
public:
	// Constructors and distructor

	Vector();
	Vector(const int n);
	Vector(const Vector& V);
	~Vector();

	// Set functions

	inline void set(const int n, const real x)
	{ _U[n-1] = x; }

	void operator= (const Vector& V);

	void resize (const int n);
	
	void raise(const int n); // raise size but reset values 

	void reset();

	// Get functions

	inline real at(const int n) const
	{ return _U[n-1]; } 

	inline int size() const
	{ return _size; }

	// Additions

	void operator+= (const Vector& V);
	Vector operator+ (const Vector& V) const;

	// Substractions

	void operator-= (const Vector& V);
	Vector operator- (const Vector& V) const;
	Vector operator- () const;

	// Multiplication with real

	void operator*= (const real x);
	
	// Division with real
	
	Vector operator/ (const real x);

	// Dot product

	real operator* (const Vector& V) const;

private:
	real* _U;
	int _size;
};

/*
_______________________________________________________________________________

External functions
_______________________________________________________________________________

*/

Vector 	operator* (const real x, const Vector& V);
real 	max(const Vector& V);

inline int dim(const Vector& V)
{ return V.size();}
