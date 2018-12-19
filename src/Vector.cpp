#include <iostream>
#include <stdlib.h>
using namespace std;
#include "Global.h"
#include "Vector.h"


/*
_______________________________________________________________________________

Empty constructor
_______________________________________________________________________________

*/

Vector::Vector()
{
	_U = new real;
	_size = 1;
	*_U = 0.;
}

/*
_______________________________________________________________________________

Constructor with size
_______________________________________________________________________________

*/

Vector::Vector(const int n)
{
	int i;

	_U = new real[n];
	_size = n;
	
	for (i = 1; i <= n; i++)
	  set(i, 0.);
	
}

/*
_______________________________________________________________________________

Copy constructor
_______________________________________________________________________________

*/

Vector::Vector(const Vector& V)
{
	int i;
	
	_U = new real[V.size()];
	_size = V.size();
  
	for (i = 1; i <= _size; i++)
	  set(i,V.at(i));
}

/*
_______________________________________________________________________________

Distructor
_______________________________________________________________________________

*/

Vector::~Vector()
{
    delete[] _U;
}

/*
_______________________________________________________________________________

Reset
_______________________________________________________________________________

*/

void Vector::reset()
{
	int i;

	for (i=1;i<=_size;i++)
	  set(i,0.);
}

/*
_______________________________________________________________________________

Resize
_______________________________________________________________________________

*/

void Vector::resize(const int n)
{
  Vector V(size());
  int minSize;
  int i;
  
  // Store local vector into V
  
  for (i = 1; i <= _size; i++)
    V.set(i,at(i));
  
  // Resize local vector and reset values
  
  raise(n);
  
  // Affect former values of the local vector
  
  minSize = min(V.size(),size());
  for (i = 1; i <= minSize; i++)
    set(i,V.at(i));
}
/*
_______________________________________________________________________________

Raise vector (with empty values)
_______________________________________________________________________________

*/

void Vector::raise(const int n)
{
  delete[] _U;
  _U = new real[n];
  _size = n;
  reset();
}

/*
_______________________________________________________________________________

Operator equality (assessment)
_______________________________________________________________________________

*/

void Vector::operator= (const Vector& V)
{
  int i;
  
  if (size() != V.size())
    raise(V.size());
  
  for (i = 1; i <= _size; i++)
    set(i,V.at(i));
  
}

/*
_______________________________________________________________________________

Add a vector into the same vector
_______________________________________________________________________________

*/

void Vector::operator+= (const Vector& V)
{
	int i;

#ifdef debug_mode
	if (size() != V.size())
	{
		cout << "in Vector.cpp: operator += \n";
		cout << "Vectors have different sizes \n";
		exit(-1);
	}
#endif

	for (i=1; i<=size(); i++)
		set(i,at(i)+V.at(i));	
}

/*
_______________________________________________________________________________

Add a vector into another vector
_______________________________________________________________________________

*/

Vector Vector::operator+ (const Vector& V) const
{
	int i;
	Vector W(size());

#ifdef debug_mode
	if (size() != V.size())
	{
		cout << "in Vector.cpp: operator + \n";
		cout << "Vectors have different sizes \n";
		exit(-1);
	}
#endif

	for (i=1; i<=size(); i++)
		W.set(i,at(i)+V.at(i));	

	return W;
}

/*
_______________________________________________________________________________

Substract a vector into the same vector
_______________________________________________________________________________

*/

void Vector::operator-= (const Vector& V)
{
	int i;

#ifdef debug_mode
	if (size() != V.size())
	{
		cout << "in Vector.cpp: operator -= \n";
		cout << "Vectors have different sizes \n";
		exit(-1);
	}
#endif

	for (i=1; i<=size(); i++)
		set(i,at(i)-V.at(i));	
}

/*
_______________________________________________________________________________

Substract a vector into another vector
_______________________________________________________________________________

*/

Vector Vector::operator- (const Vector& V) const
{
	int i;
	Vector W(size());

#ifdef debug_mode
	if (size() != V.size())
	{
		cout << "in Vector.cpp: operator - \n";
		cout << "Vectors have different sizes \n";
		exit(-1);
	}
#endif

	for (i=1; i<=size(); i++)
		W.set(i,at(i)-V.at(i));	

	return W;
}

/*
_______________________________________________________________________________

Opposite of a vector
_______________________________________________________________________________

*/

Vector Vector::operator- () const
{
  int i;
  Vector W(size());
  
  for (i=1; i<=size(); i++)
    W.set(i,-at(i));	

  return W;
 
}

/*
_______________________________________________________________________________

Multiply a vector with a real
_______________________________________________________________________________

*/

void Vector::operator*= (const real x)
{
	int i;

	for (i=1; i<=size(); i++)
		set(i,at(i)*x);	
}

/*
_______________________________________________________________________________

Divide a vector with a real
_______________________________________________________________________________

*/

Vector Vector::operator/ (const real x)
{
	int i;
	Vector result(size());

	for (i=1; i<=size(); i++)
		result.set(i,at(i)/x);	

	return result;
}

/*
_______________________________________________________________________________

Dot product
_______________________________________________________________________________

*/

real Vector::operator* (const Vector& V) const
{
	int i;
	real result=0.;

	
#ifdef debug_mode
	if (size() != V.size())
	{
		cout << "in Vector.cpp: operator* (dot product) \n";
		cout << "Vectors have different sizes \n";
		exit(-1);
	}
#endif

	for (i=1; i<=size(); i++)
		result += at(i)*V.at(i);	

	return result;
}

/*
_______________________________________________________________________________

EXTERNAL: Multiply a vector with a real (store into a new real)
_______________________________________________________________________________

*/

Vector operator* (const real x, const Vector& V)
{
	int i;
	Vector W(V.size());

	for (i=1; i<=V.size(); i++)
		W.set(i,V.at(i)*x);

	return W;

}

/*
_______________________________________________________________________________

EXTERNAL: Maximum of a vector
_______________________________________________________________________________

*/

real max(const Vector& V)
{
    int i;
    real result=0.;
    
    for (i=1; i<=V.size();i++)
      result=MAX(result,ABS(V.at(i)));
    
    return result;
}